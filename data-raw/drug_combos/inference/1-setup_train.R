library(crossmeta)
library(ccdata)
library(RcppCNPy)

source('~/Documents/Batcave/GEO/ccdata/data-raw/drug_combos/combo_utils.R')


# Get Input Genes (combo data) ------------------------------

setwd("~/Documents/Batcave/GEO/ccdata/data-raw/drug_combos/inference")
data_dir <- paste("~/Documents/Batcave/GEO/1-meta", "COMBOS", sep="/")

# load previous analyses
gse_names <- read.table("gse_names.csv", quote="\"",
                        comment.char="", stringsAsFactors=FALSE)$V1

anals <- load_diff(gse_names, data_dir, "SYMBOL")

# get all genes
genes <- lapply(anals, function(anal) {
    row.names(anal$top_tables[[1]])
})

# get input (common) genes
input_genes <- Reduce(intersect, genes)
saveRDS(input_genes, "input_genes.rds")


# Add Training Data (from '1-meta') ------

# get GSEs with most HGU133A genes
data(cmap_es)
hgu_genes <- row.names(cmap_es)

all_paths <- list.files("~/Documents/Batcave/GEO/1-meta",
                        "_diff_expr_symbol.rds",
                        recursive = TRUE, full.names = TRUE)

train_paths <- c()
genes <- list()
num_cons  <- 0   # number of contrasts in train_paths esets

for (i in seq_along(all_paths)) {
    eset <- readRDS(all_paths[i])
    miss <- setdiff(hgu_genes, row.names(eset$top_tables[[1]]))

    if (length(miss) <= 1000) {
        j <- length(genes) + 1

        train_paths <- c(train_paths, all_paths[i])
        genes[[j]] <- row.names(eset$top_tables[[1]])
        num_cons <- num_cons + length(eset$top_tables)
    }
    print(i)
}

# number of common genes in train_paths esets
all_genes <- Reduce(intersect, genes)

saveRDS(all_genes, "all_genes.rds")
saveRDS(train_paths,  "train_paths.rds")

# update input genes to ensure in all genes
input_genes <- readRDS("input_genes.rds")
input_genes <- intersect(all_genes, input_genes)

saveRDS(input_genes, "input_genes.rds")


# Setup Training Data --------------------

data(cmap_es)
train_paths <- readRDS("train_paths.rds")
input_genes <- readRDS("input_genes.rds")
all_genes   <- readRDS("all_genes.rds")

# load analyses and add dprimes to top tables
anals <- lapply(train_paths, readRDS)
anals <- add_es(anals, "dprime")

# put all top tables into a list
train <- unlist(lapply(anals, `[[`, 'top_tables'), recursive = FALSE)

# remove non-target/input genes and non-dprime columns
train <- lapply(train, `[`, all_genes, 'dprime')
train <- as.data.frame(train, row.names = all_genes)

# combine with cmap_es
train <- cbind(cmap_es[all_genes, ], train)

# seperate into X and y
X <- train[input_genes, ]
y <- train[setdiff(all_genes, input_genes), ]

# save
X <- t(as.matrix(X))
y <- t(as.matrix(y))

saveRDS(X, "data/X.rds")
saveRDS(y, "data/y.rds")
npySave('data/X.npy', X, mode='w')
npySave('data/y.npy', y, mode='w')
