library(crossmeta)
library(data.table)
library(RcppCNPy)

setwd("~/Documents/Batcave/GEO/ccdata/data-raw/drug_combos")
source('combo_utils.R')


# Get Test Data (combo data) ----------------

data_dir <- paste("~/Documents/Batcave/GEO/1-meta", "COMBOS", sep="/")

input_genes <- readRDS("inference/input_genes.rds")
all_genes   <- readRDS("inference/all_genes.rds")


# load previous analyses
gse_names <- read.table("gse_names.csv", quote="\"",
                        comment.char="", stringsAsFactors=FALSE)$V1

anals <- load_diff(gse_names, data_dir, "SYMBOL")
anals <- add_es(anals)

# get top tables
test <- lapply(anals, `[[`, 'top_tables')
test <- unlist(test, recursive = FALSE)


# seperately keep track of vardprime values
vars <- lapply(test, function (x) {
    x[row.names(x) %in% all_genes, 'vardprime', drop = FALSE]
})

# keep only target/input (all) genes and dprime values
test <- lapply(test, function (x) {
    x[row.names(x) %in% all_genes, 'dprime', drop = FALSE]
})


# Setup Test Data (combo data) ----------------

# merge all test data
all_vars  <- lapply(vars, as.data.table, keep.rownames = TRUE)
all_data  <- lapply(test, as.data.table, keep.rownames = TRUE)
all_genes <- data.table(rn = all_genes)

all_data <- lapply(all_data, function(x) {
    merge(all_genes, x, all.x = TRUE, on = "rn")
})

all_vars <- lapply(all_vars, function(x) {
    merge(all_genes, x, all.x = TRUE, on = "rn")
})

test <- lapply(all_data, function(x) {
    x$dprime
})

vars <- lapply(all_vars, function(x) {
    x$vardprime
})

test <- as.data.frame(test)
vars <- as.data.frame(vars)
test <- t(test)
vars <- t(vars)
colnames(test) <- all_data[[1]]$rn
colnames(vars) <- all_vars[[1]]$rn

# determine complete samples (don't need preds)
cplt <- apply(test, 1, function(x) sum(is.na(x)) == 0)

# get order of genes in X and y
X <- readRDS("inference/data/X.rds")
y <- readRDS("inference/data/y.rds")

Xtest <- test[!cplt, colnames(X)]
ytest <- test[!cplt, colnames(y)]

Xvars <- vars[!cplt, colnames(X)]
yvars <- vars[!cplt, colnames(y)]

# will need complete samples to train for combo predictions
cpltd <- test[cplt, c(colnames(X), colnames(y))]
cpltv <- vars[cplt, c(colnames(X), colnames(y))]

# save Xtest, ytest, and cplt
saveRDS(Xtest, "inference/data/Xtest.rds")
saveRDS(ytest, "inference/data/ytest.rds")
saveRDS(Xvars, "inference/data/Xvars.rds")
saveRDS(yvars, "inference/data/yvars.rds")
saveRDS(cpltd, "inference/data/cpltd.rds")
saveRDS(cpltv, "inference/data/cpltv.rds")

npySave('inference/data/Xtest.npy', Xtest, mode='w')
npySave('inference/data/ytest.npy', ytest, mode='w')
