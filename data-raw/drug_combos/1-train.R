library(crossmeta)
library(ccdata)

setwd("~/Documents/Batcave/GEO/ccdata/data-raw/drug_combos")

# Load Data -------------------------

data_dir <- paste("~/Documents/Batcave/GEO/1-meta", "COMBOS", sep="/")

# load previous analyses
gse_names <- read.table("gse_names.csv", quote="\"",
                        comment.char="", stringsAsFactors=FALSE)$V1

anals <- load_diff(gse_names, data_dir, "SYMBOL")

#load esets
esets <- load_raw(gse_names, data_dir)


#GSE17672.GPL3921 incomplete
#esets$"GSE17672.GPL3921" <- NULL


# Differential Expression -------------------------

#reload previous analysis
prev <- load_diff(gse_names, data_dir, "SYMBOL")

#re-run previous analysis
anals <- diff_expr(esets, data_dir, "PROBE")


# Setup Training Data -------------------------


combo_train <- readRDS("combo_train.rds")
selections  <- readRDS("selections.rds")
combo_train <- setup_combo_data(anals, selections)

saveRDS(combo_train$data, "combo_train.rds")
saveRDS(combo_train$selections, "selections.rds")
