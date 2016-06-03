library(xgboost)

setwd("/home/alex/Documents/Batcave/GEO/ccdata/data-raw/drug_combos")
# Setup -------------------------

#load data
train <- readRDS("combo_train.rds")

#get label
y <- train$combo_dprime

#remove combo info from X
cb_cols <- grepl("combo", colnames(train))
ls_cols <- grepl("dprime|adj.P.Val", colnames(train))
X <- train[, ls_cols & !cb_cols]

#swap values between drug1 and drug2 columns
d1_cols <- grep("drug1", colnames(X), value=T)
d2_cols <- grep("drug2", colnames(X), value=T)

Xr <- X
colnames(Xr) <- c(d2_cols, d1_cols)
Xr <- Xr[, colnames(X)]

#for each row, randomly use swapped or original orientation
filt <- sample(c(TRUE, FALSE), nrow(X), replace=TRUE)
X[filt, ] <- Xr[filt, ]


# Grid Search -------------------------

evalerror <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    err <- 1 - (sum(sign(preds) == sign(labels)) / length(labels))
    return(list(metric = "error", value = err))
}


grid <- expand.grid(subsample=c(1, 0.75, 0.5),
                    colsample_bytree=1,
                    max.depth=c(10, 8, 6, 4))
ntrees <- 1000

#Build a xgb.DMatrix object
dtrain <- xgb.DMatrix(data=as.matrix(X), label=y)

#cleanup
rm(train, Xr, filt, d1_cols, d2_cols, cb_cols, X, y)
gc()

res <- data.frame(test_error=numeric(0),
                  train_error=numeric(0),
                  best_iter=numeric(0))

for (i in 1:nrow(grid)) {

    #Extract Parameters to test
    sub <- grid[["subsample"]][i]
    col <- grid[["colsample_bytree"]][i]
    dep <- grid[["max.depth"]][i]

      params <- list(objective           = "reg:linear",
                     eta                 = 0.2,
                     subsample           = 1,
                     colsample_bytree    = 1,
                     max_depth           = 8,
                     tree_method         = "approx",
                     eval_metric         = evalerror)

      history <- xgb.cv(params, dtrain, 580,
                        nfold = 5, early.stop.round = 50,
                        maximize = FALSE, prediction = TRUE)


    #Save rmse/best iteration
    test_error  <- min(history$test.error.mean)
    train_error <- min(history$train.error.mean)
    best_iter   <- which.min(history$test.error.mean)

    params <- paste(sub, col, dep, sep="-")
    res[params, ] <- c(test_error, train_error, best_iter)

    cat("\n")
    print(res)
    cat("\n")

    rm(history)
    gc()
}

#            test_error  train_error   best_iter
# 1-1-10      0.233265    0.229447       398
# 0.75-1-10   0.234184    0.232794       951
# 0.5-1-10    0.234597    0.233849       484
# 1-1-8       0.233480    0.231898       580   <- THIS
# 0.75-1-8    0.234885    0.234622       395
# 0.5-1-8     0.234982    0.234721       442
# 1-1-6       0.234197    0.233761       649
# 0.75-1-6    0.235096    0.235030       457
# 0.5-1-6     0.235246    0.235037       555
# 1-1-4       0.234859    0.234727       916
# 0.75-1-4    0.235377    0.235246       822
# 0.5-1-4     0.235484    0.235394       871





# Model -------------------------

combo_model <- xgboost(data=dtrain, nround=580,
                       objective = "reg:linear", eta=0.2,
                       subsample=1, colsample_bytree=1, max.depth=8, tree_method='exact')

#plot feature importance
imp_matrix <- xgb.importance(feature_names=colnames(X), model=combo_model)
print(xgb.plot.importance(importance_matrix=imp_matrix))

#save model
devtools::use_data(combo_model, overwrite = TRUE)



# Analyse Model -------------------

#preds data.frame
pdf <- data.frame(xgb = history$pred,
                  scaled = history$pred * 1.623,
                  avg = (train$drug1_dprime + train$drug2_dprime) / 2,
                  actual = y)

#line of best fit for actual vs xgb (and scaled preds)
fitpr <- lm(xgb ~ actual , data=pdf)
fitsc <- lm(scaled ~ actual , data=pdf)


#plot actual vs predicted
library(ggplot2)

ggplot(pdf) +
    #scaled vs actual hexbins
    geom_hex(aes(x = scaled, y = actual)) +
    #xgb vs actual line
    geom_abline(intercept = coef(fitpr)[1], slope = coef(fitpr)[2], colour = "red") +
    #scaled vs actual line
    geom_abline(intercept = coef(fitsc)[1], slope = coef(fitsc)[2], colour = "black") +
    #perfect line (x = y)
    geom_abline(intercept = 0, slope = 1, colour = "green")


#plot residuals
pdf$resid_pr <- pdf$actual - pdf$xgb
pdf$resid_sc <- pdf$actual - pdf$scaled

ggplot(pdf) + geom_hex(aes(x = xgb, y = resid_pr))
ggplot(pdf) + geom_hex(aes(x = scaled, y = resid_sc))

#plot prediction accuracy as a function of absolute effect size
cuts <- seq(0, 6.1, 0.1)  # 61 bins

get_accuracy <- function(preds, y, cuts) {

    acs <- c()
    for (i in 1:61) {

        #define bin range
        bin <- abs(y) > cuts[i] & abs(y) <= cuts[i+1]

        #get accuracy in bin
        ac <- sum(sign(y[bin]) == sign(preds[bin])) / sum(bin)
        acs <- c(acs, ac)
    }
    return(acs)
}

xgb_accuracy <- get_accuracy(pdf$xgb, y, cuts)
avg_accuracy <- get_accuracy(pdf$avg, y, cuts)


adf <- data.frame(accuracy = c(xgb_accuracy, avg_accuracy),
                  Model = c(rep(c("XGBoost", "Average"), each = length(avg_accuracy))),
                  es = rep(cuts[-62], 2))


ggplot(adf) +
    geom_line(aes(x = es, y = accuracy, colour = Model)) +
    ylab("Classification Accuracy") +
    xlab("Absolute Effect Size") +
    scale_x_continuous(breaks = 0:6) +
    scale_y_continuous(breaks = seq(0.5, 1, 0.05))

#overall accuracy
sum(sign(y) == sign(pdf$xgb)) / length(y)

# 0.769
# 0.760 for 'avg' model
# 0.708 for 'drug1 or drug2' model
