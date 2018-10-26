library(crossmeta)
library(ccdata)
library(preprocessCore)
library(RcppCNPy)

setwd("~/Documents/Batcave/GEO/ccdata/data-raw/drug_combos/inference")
data_dir <- paste("~/Documents/Batcave/GEO/1-meta", "COMBOS", sep="/")

# load previous analyses
gse_names <- read.table("gse_names.csv", quote="\"",
                        comment.char="", stringsAsFactors=FALSE)$V1


# Adjust for SVs ------------------------------

# diff_anal in crossmeta: exprs(eset) <- exprs_sva

esets <- load_raw(gse_names, data_dir)
prev <- load_diff(gse_names, data_dir, "SYMBOL")
anals <- diff_expr(esets, data_dir, "SYMBOL", prev)

# Keep HGU133A genes ---------------------------

anals <- load_diff(gse_names, data_dir, "SYMBOL")

# remove genes not in cmap_es
dat <- lapply(anals, function(x) exprs(x$eset))
data(cmap_es)
genes <- row.names(cmap_es)

dat <- lapply(dat, function(x) {
    x[row.names(x) %in% genes, ]
})


# Quantile Normalize HGU133A genes --------------

# load 10,000 random samples from cleaned training data
set.seed(0)
ids <- sample(113927, 10000)

y <- readRDS("clean/y.rds")[ids, ] ; gc()
X <- readRDS("clean/X.rds")[ids, ] ; gc()

train <- cbind(X, y) ; gc()
train <- t(train) ; gc()     # genes in rows (same as dat)

# quantile normalize each combo dataset to train with same genes (~30min)
dat <- lapply(dat, function(x) {
    target <- as.vector(train[row.names(train) %in% row.names(x), ])
    normalize.quantiles.use.target(x, target, copy = FALSE)
})

saveRDS(dat, "qnorm_combos.rds")


# Setup Test Data --------------

dat <- readRDS("qnorm_combos.rds")

# get X genes in right order
X_genes <- colnames(readRDS("clean/X.rds")) ; gc()

# subset, order, transpose, scale, and save each matrix
Xt <- lapply(dat, function(mt) {
    t(mt[X_genes, ])
})

Xt <- Reduce(rbind, Xt)

Xt <- scale(Xt)


