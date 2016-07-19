library(crossmeta)
library(ccdata)

setwd("~/Documents/Batcave/GEO/ccdata/data-raw/drug_combos")

# Load Data -------------------------

data_dir <- paste("~/Documents/Batcave/GEO/1-meta", "COMBOS", sep="/")

#load gse names from csv

get_raw(gse_names, data_dir)

#load esets
esets <- load_raw(gse_names, data_dir)

#backup
esets_copy <- esets

#GSE17672.GPL3921 incomplete
#esets$"GSE17672.GPL3921" <- NULL


# Differential Expression -------------------------

#reload previous analysis
anals <- load_diff(gse_names, data_dir, "PROBE")

#re-run previous analysis
anals <- diff_expr(esets, data_dir, "PROBE", prev)


# Setup Training Data -------------------------


combo_train <- readRDS("combo_train.rds")
selections  <- readRDS("selections.rds")
combo_train <- setup_combo_data(anals, selections)

saveRDS(combo_train$data, "combo_train.rds")
saveRDS(combo_train$selections, "selections.rds")
