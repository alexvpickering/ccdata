library(xgboost)
library(ccmap)
library(ccdata)


# Setup -------------------------

#load data
data(combo_train)
train <- combo_train$data

#get label (1 if up)
y <- ifelse(train$combo_dprime > 0, 1, 0)

#remove combo info from X
cb_cols <- grepl("combo", colnames(train))
X <- train[, !cb_cols]

#swap values between drug1 and drug2 columns
d1_cols <- grep("drug1", colnames(X), value=T)
d2_cols <- grep("drug2", colnames(X), value=T)

Xr <- X
colnames(Xr) <- c(d2_cols, d1_cols)
Xr <- Xr[, colnames(X)]

#for each row, randomly use swapped or original orientation
filt <- sample(c(TRUE, FALSE), nrow(X), replace=TRUE)
X[filt, ] <- Xr[filt, ]


#cleanup
rm(train, Xr, filt, d1_cols, d2_cols, cb_cols)


# Grid Search -------------------------

grid <- expand.grid(subsample=c(1, 0.75, 0.5), 
                     colsample_bytree=c(1, 0.75, 0.5),
                     max.depth=c(10, 7, 4))
ntrees <- 2000

#Build a xgb.DMatrix object
dtrain <- xgb.DMatrix(data=as.matrix(X), label=y)

hyps <- c()
for (i in 1:nrow(grid)) {
 
   #Extract Parameters to test
   sub <- grid[["subsample"]][i]
   col <- grid[["colsample_bytree"]][i]
   dep <- grid[["max.depth"]][i]
 
   history <- xgb.cv(data=dtrain, nrounds=ntrees, nfold=5, 
                   objective="binary:logistic", eta=0.1,                               
                   subsample=sub, colsample_bytree=col, max.depth=dep,
                   early.stop.round=50, maximize=F)
 
   history <- as.data.frame(history)
 
   #Save error of the best iteration
   error <- min(history$test.error.mean)
   names(error) <- paste(sub, col, dep, sep="-")
 
   cat ("\nerror:", error, "sub-col-dep:", sub, col, dep, "\n\n")
   hyps <- c(hyps, error)
 }

# error: 0.218157 
# sub-col-dep: 0.75 1 10 
# best iteration: 480


# Model -------------------------

combo_model <- xgboost(data=dtrain, nround=480, 
                       objective = "binary:logistic", eta=0.1, 
                       subsample=0.75, colsample_bytree=1, max.depth=10)

#plot feature importance
imp_matrix <- xgb.importance(feature_names=colnames(X), model=combo_model)
print(xgb.plot.importance(importance_matrix=imp_matrix))

#save model
devtools::use_data(combo_model)