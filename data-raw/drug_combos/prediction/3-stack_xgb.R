library(xgboost)
library(data.table)

setwd("/home/alex/Documents/Batcave/GEO/ccdata/data-raw/drug_combos/prediction")

# load data to train stacker on
train <- fread("data/stack/train.csv")
names(train) <- c("net_preds", "drug1_dprime", "drug2_dprime", "drug1_vardprime",
                  "drug2_vardprime", "combo_dprime")

X <- train[, !"combo_dprime", with=FALSE]

y <- train$combo_dprime

# evaluation metric
accuracy <- function(preds, dtrain) {
    labels <- getinfo(dtrain, "label")
    acc <- sum(sign(preds) == sign(labels)) / length(labels)
    return(list(metric = "acc", value = acc))
}


# grid <- expand.grid(subsample=c(1, 0.75, 0.5),
#                     colsample_bytree=c(1, 0.75, 0.5),
#                     max.depth=c(10, 8, 6, 4))

#Build a xgb.DMatrix object
dtrain <- xgb.DMatrix(data=as.matrix(X), label=y, missing = NA)

# res <- data.frame(test_error=numeric(0),
#                   train_error=numeric(0),
#                   best_iter=numeric(0))
#
# for (i in 1:nrow(grid)) {
#
#     #Extract Parameters to test
#     sub <- grid[["subsample"]][i]
#     col <- grid[["colsample_bytree"]][i]
#     dep <- grid[["max.depth"]][i]
#
#     params <- list(objective           = "reg:linear",
#                    eta                 = 0.2,
#                    subsample           = sub,
#                    colsample_bytree    = col,
#                    max_depth           = dep,
#                    eval_metric         = evalerror)
#
#     history <- xgb.cv(params, dtrain, 600,
#                       nfold = 5, early.stop.round = 20,
#                       maximize = FALSE, prediction = FALSE)
#
#
#     #Save rmse/best iteration
#     test_error  <- min(history$test.error.mean)
#     train_error <- min(history$train.error.mean)
#     best_iter   <- which.min(history$test.error.mean)
#
#     params <- paste(sub, col, dep, sep="-")
#     res[params, ] <- c(test_error, train_error, best_iter)
#
#     cat("\n")
#     print(res)
#     cat("\n")
#
#     rm(history)
#     gc()
# }


#               test_error   train_error    best_iter     included_cols
#    1-1-15     0.198213     0.192141         8           all
#    1-1-10     0.197294     0.193615       578           all
# 0.75-1-10     0.208687     0.206179       592           all - net_preds
# 0.75-1-10     0.202828     0.202300        24           all - vardprime
#    1-1-10     0.197571     0.194837       498           all - average
# 0.75-1-10     0.210437     0.210042        25           dprimes

# avg   : 0.210365 (2424 genes)
# +vars : 0.208687 (- 19 genes)
# +nets : 0.202828 (- 86 genes)
# +both : 0.198203 (-140 genes)

# Model -------------------------

param_updates <- list(eta = c(0.5, 0.5, rep(0.15, 6)))

history <- xgb.cv(data = dtrain, nround = 8,
                  objective = "reg:linear", eta = 0.5,
                  max.depth = 15, nfold = 5,
                  prediction = TRUE, feval = accuracy, callbacks = param_updates)

saveRDS(history, "data/history.rds")

xgb_mod <- xgboost(data=dtrain, nround=8,
                       objective = "reg:linear", eta=0.5,
                       max.depth=15, callbacks = param_updates)

devtools::use_data(xgb_mod, overwrite = TRUE)

#plot feature importance
imp_matrix <- xgb.importance(feature_names=colnames(X), model=xgb_mod)
print(xgb.plot.importance(importance_matrix=imp_matrix))
