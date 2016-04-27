library(crossmeta)
library(ccdata)

# Load Data -------------------------

#GSE2487, GSE2450, GSE30644 - - > singular

data_dir <- paste("~/Documents/Batcave/GEO/1-meta", "data", "COMBOS", sep="/")

affy_names <- c("GSE49416", "GSE12791", "GSE44857", "GSE38972", "GSE35218", "GSE30524",
                "GSE5681", "GSE26599", "GSE28646", "GSE22445", "GSE73906", "GSE55624",
                "GSE66854", "GSE66837", "GSE35696", "GSE7500", "GSE57915", "GSE54169",
                "GSE55541", "GSE51717", "GSE20405", "GSE43685", "GSE46362", "GSE48258",
                "GSE36572", "GSE38663", "GSE31683", "GSE22025", "GSE35230", "GSE33562",
                "GSE33366", "GSE16179", "GSE26114", "GSE28896", "GSE27973", "GSE22589",
                "GSE20876", "GSE20115", "GSE9649", "GSE8615", "GSE11552", "GSE7568",
                "GSE5054", "GSE7035", "GSE6960", "GSE6914", "GSE72359", "GSE9988",
                "GSE14964", "GSE17672", "GSE28448", "GSE71347")


get_raw(affy_names, data_dir)

#load esets
esets <- load_raw(affy_names, data_dir)

#backup
esets_copy <- esets

#GSE17672.GPL3921 incomplete
esets$"GSE17672.GPL3921" <- NULL

esets <- commonize(esets, "PROBE")



# Differential Expression -------------------------


#run initial analysis
anals <- diff_expr(esets, data_dir, "PROBE")

#reload previous analysis
anals <- load_diff(affy_names, data_dir, probe=TRUE)

#backup
anals_copy <- anals

#re-run previous analysis
anals <- diff_expr(esets, data_dir, "PROBE", anals)


# Setup Training Data -------------------------


#not individual & combo: GSE12791, GSE35218
diff_exprs$GSE12791 <- NULL
diff_exprs$GSE35218 <- NULL

data(combo_train)

combo_train <- setup_combo_data(diff_exprs, combo_train$selections)

devtools::use_data(combo_train)
