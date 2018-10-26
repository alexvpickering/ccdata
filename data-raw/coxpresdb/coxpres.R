setwd("~/Documents/Batcave/GEO/ccdata/data-raw/coxpresdb")


data_dir <- "~/Documents/Batcave/GEO/ccdata/data-raw/coxpresdb/Hsa/Hsa.v13-01.G20280-S73083.rma.mrgeo.d"
files <- list.files(data_dir)
paths <- file.path(data_dir, files)

library(data.table)
hs_coxpres <- data.table(gene1 = rep(as.character(files), each = 200),
                         gene2 = "",
                         mrank = 0)


data <- lapply(paths, fread, skip = 1, nrows = 200)
data <- rbindlist(data)

hs_coxpres$gene2 <- as.character(data$V1)
hs_coxpres$mrank <- data$V2

saveRDS(hs_coxpres, "hs_coxpres.rds")


