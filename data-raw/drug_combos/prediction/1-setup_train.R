library(crossmeta)
library(ccdata)
library(RcppCNPy)
library(data.table)

setwd("~/Documents/Batcave/GEO/ccdata/data-raw/drug_combos")
source('combo_utils.R')

# Load Data -------------------------

data_dir <- paste("~/Documents/Batcave/GEO/1-meta", "COMBOS", sep="/")

# load previous analyses
gse_names <- read.table("gse_names.csv", quote="\"",
                        comment.char="", stringsAsFactors=FALSE)$V1

anals <- load_diff(gse_names, data_dir, "SYMBOL")

# get selections for analyses
selections <- readRDS("selections.rds")
#selections <- select_combo_data(anals, selections)

# load predicted values
Xtest <- readRDS("inference/data/Xtest.rds")
ytest <- readRDS("inference/data/ytest.rds")
preds <- read.table("inference/data/preds.csv", sep=",")

Xvars <- readRDS("inference/data/Xvars.rds")
yvars <- readRDS("inference/data/yvars.rds")
vars  <- cbind(Xvars, yvars)

colnames(preds) <- c(colnames(Xtest), colnames(ytest))
row.names(preds) <- row.names(Xtest)

# load complete cases and combine
cpltd <- readRDS("inference/data/cpltd.rds")
cpltv <- readRDS("inference/data/cpltv.rds")

cbo_data <- rbind(preds, cpltd)
cbo_vars <- rbind(vars,  cpltv)

row.names(cbo_data) <- gsub("^GSE\\d+.GSE|^GSE\\d+.GPL\\d+.GSE",
                            "GSE", row.names(cbo_data))

row.names(cbo_vars) <- gsub("^GSE\\d+.GSE|^GSE\\d+.GPL\\d+.GSE",
                            "GSE", row.names(cbo_vars))

# setup train data
train  <- combine_combo_data(cbo_data, selections)
trainv <- combine_combo_data(cbo_vars, selections)

X <- train$X
y <- train$y
Xv <- trainv$X
yv <- trainv$y

# for X, swap values between drug1 and drug2 columns
d1_cols <- grep("drug1", colnames(X), value=TRUE)
d2_cols <- grep("drug2", colnames(X), value=TRUE)

Xr  <- X
Xvr <- Xv

colnames(Xr)  <- c(d2_cols, d1_cols)
colnames(Xvr) <- c(d2_cols, d1_cols)

Xr  <- Xr[, colnames(X)]
Xvr <- Xvr[, colnames(X)]

# for each row of X, randomly use swapped or original orientation
set.seed(0)
filt <- sample(c(TRUE, FALSE), nrow(X), replace=TRUE)
X[filt, ]  <- Xr[filt, ]
Xv[filt, ] <- Xvr[filt, ]

saveRDS(X, "prediction/data/X.rds")
saveRDS(y, "prediction/data/y.rds")
saveRDS(Xv, "prediction/data/Xv.rds")
saveRDS(yv, "prediction/data/yv.rds")

npySave('prediction/data/X.npy', X, mode='w')
npySave('prediction/data/y.npy', y, mode='w')
npySave('prediction/data/Xv.npy', Xv, mode='w')
npySave('prediction/data/yv.npy', yv, mode='w')

# get rank of absolute gene expression for each sample within y
y <- t(apply(abs(y), 1, frank, ties.method = 'dense'))

saveRDS(y, "prediction/data/y_ranked.rds")
npySave('prediction/data/y_ranked.npy', y, mode='w')

