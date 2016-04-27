library(reshape2)

setwd("~/Documents/Batcave/GEO/ccdata/data-raw/atc4")

#load in cmap drugs
cmap_instances <- read.table("~/Documents/Batcave/GEO/ccdata/data-raw/raw/cmap_instances_02.csv",
                             header=T, sep="\t", quote='',
                             fill=T, stringsAsFactors=F)

drugs <- unique(cmap_instances$cmap_name)

#load in atc codes
load("AllData-WHOCC-dump-2016-02-12.RData")
atc <- as.data.frame(AllData[["atc"]], stringsAsFactors=F)

#get atc codes for drugs (7 chars) in cmap
atc <- atc[nchar(atc$key) == 7, ]
atc <- atc[atc$name %in% drugs, ]

#obtain list of 4th level atc codes
atc4 <- list()
for (drug in unique(atc$name)) {
    #full atc codes
    keys <- atc[atc$name == drug, ]$key
    #4th level atc codes
    keys <-gsub("(.+)\\d\\d", "\\1", keys)
    atc4[[drug]] <- unique(keys)
}

devtools::use_data(atc4, ccdata)
