library(crossmeta)
library(data.table)
library(hgu133a.db)
library(RcppCNPy)

setwd("~/Documents/Batcave/GEO/ccdata/data-raw/drug_combos/inference")
data_dir <- paste("~/Documents/Batcave/GEO/1-meta", "COMBOS", sep="/")

source("io.R")

# Common Genes (combo data) ------------------------------

# load previous analyses
gse_names <- read.table("gse_names.csv", quote="\"",
                        comment.char="", stringsAsFactors=FALSE)$V1

anals <- load_diff(gse_names, data_dir, "SYMBOL")

# get all genes
all_genes <- lapply(anals, function(anal) {
    row.names(anal$top_tables[[1]])
})

# get common genes
common_genes <- Reduce(intersect, all_genes)
saveRDS(common_genes, "common_genes.rds")

# Unique Samples (columns of bgedv2) ------------------------------

# load subset of probes
set.seed(0)
rows <- sample(22268, 6000)

ds <- parse.gctx("raw/bgedv2_QNORM.gctx", rid = rows)@mat
ds <- t(ds); gc()
ds <- as.data.table(ds); gc()
dups <- duplicated(ds)

unique_samples <- which(!dups)
saveRDS(unique_samples, "unique_samples.rds")

rm(ds); gc()

# Probe IQRs (rows of bgedv2) ------------------------------
cols <- readRDS("unique_samples.rds")

x <- 1:22268
n <- 10
bins <- split(x, sort(x%%n))

iqrs <- c()

for (bin in bins){
    ds <- parse.gctx("raw/bgedv2_QNORM.gctx", rid = bin, cid = cols)@mat
    iqrs <- c(iqrs, matrixStats::rowIQRs(ds))
    rm(ds); gc()
}

saveRDS(iqrs, "iqrs.rds")

# Probes to Symbols (bgedv2) ------------------------------

dt         <- data.table(iqr = readRDS("iqrs.rds"))
dt$PROBEID <- row.names(parse.gctx("raw/bgedv2_QNORM.gctx", cid = 1)@mat)
dt$row     <- 1:nrow(dt)

# map from probes to symbol
map <- AnnotationDbi::select(hgu133a.db, dt$PROBEID, "SYMBOL")
map <- as.data.table(map)
map <- map[!is.na(SYMBOL), ]
map$SYMBOL <- toupper(map$SYMBOL)

# expand 1:many mappings
dt <- merge(dt, map, all.y = TRUE, by = "PROBEID", sort = FALSE)

# for rows with duplicated symbol, keep highest IQR
dt <- dt[, .SD[which.max(iqr)], by = "SYMBOL"]
saveRDS(dt, "ps_map.rds")

# Get X and y (bgedv2) ------------------------------

cols    <- readRDS("unique_samples.rds")
ps_map  <- readRDS("ps_map.rds")
X_genes <- readRDS("common_genes.rds")

Xmap <- ps_map[SYMBOL %in% X_genes]
ymap <- ps_map[!SYMBOL %in% X_genes]

#  X:
# ---

X <- parse.gctx("raw/bgedv2_QNORM.gctx", rid = unique(Xmap$row))@mat
X <- X[, cols] ; gc()               # remove non-unique samples
X <- X[Xmap$PROBEID, ] ; gc()       # expand 1-to-many probe-to-gene
row.names(X) <- Xmap$SYMBOL ; gc()  # use symbol for row names
colnames(X)  <- NULL

# samples in rows / genes in columns
X <- t(X)

saveRDS(X, "clean/X.rds")
npySave('clean/X.npy', X, mode='w')

# scale s.t. for each gene mean = 0, sd = 1
X <- scale(X)

saveRDS(X, "clean/X_scaled.rds")
npySave('clean/X_scaled.npy', X, mode='w')

rm(X); gc()


#  y:
# ---

y <- parse.gctx("raw/bgedv2_QNORM.gctx", rid = unique(ymap$row))@mat
y <- y[, cols] ; gc()               # remove non-unique samples
y <- y[ymap$PROBEID, ] ; gc()       # expand 1-to-many probe-to-gene
row.names(y) <- ymap$SYMBOL ; gc()  # use symbol for row names
colnames(y)  <- NULL ; gc()

# samples in rows / genes in columns
y <- t(y)

saveRDS(y, "clean/y.rds")
npySave('clean/y.npy', y, mode='w')

# Rank order X and y (bgedv2) ------------------------------

# for each observation, get rank of each gene
X <- readRDS("clean/X.rds")
y <- readRDS("clean/y.rds")
train <- cbind(X, y)
rm(X, y) ; gc()

for (i in 1:nrow(train)) {
    train[i, ] <- frank(train[i, ], ties.method = 'dense')
    if (i %% 1000 == 0) print(i)
}

# seperate out into X and y matrices
X <- train[, 1:3486]
saveRDS(X, "clean/X_ranked.rds")
npySave('clean/X_ranked.npy', X, mode='w')
rm(X) ; gc()

y <- train[, -(1:3486)]
saveRDS(y, "clean/y_ranked.rds")
npySave('clean/y_ranked.npy', y, mode='w')
rm(y) ; gc()
