library(dplyr)

setwd("~/Documents/Batcave/GEO/2-cmap/data/processed/prl")

#from https://github.com/francescojm/iNRG_cMap (DRUG_PRLs.ro)
load("iorio_prls.ro")


# Setup ----------------------------


#reference for order of probe names
probe_names <- iorio_prls[, 1]

#replace probe name with position of probe
for (i in seq_along(colnames(iorio_prls))) { 
  iorio_prls[,i] <- match(probe_names, iorio_prls[,i])
}

#set row names to probe names
row.names(iorio_prls) <- probe_names




# Annotate --------------------------


#map from probe names to SYMBOL
library(hgu133a.db)
map <- AnnotationDbi::select(hgu133a.db, probe_names, "SYMBOL")
map <- map[!is.na(map$SYMBOL), ]
iorio_prls <- iorio_prls[map$PROBEID, ] 

#where SYMBOL duplicated, keep SYMBOL with least central rank
cmap_prls <- data.frame(SYMBOL=unique(map$SYMBOL), stringsAsFactors=F)
middle <- nrow(iorio_prls) / 2

drug_names <- names(iorio_prls)
iorio_prls[,"SYMBOL"] <- map$SYMBOL

for (drug in drug_names) {
  #get drug prl, SYMBOL, and distance from center
  cmap_prl <- iorio_prls[, c(drug, "SYMBOL")]
  cmap_prl$dist <- abs(cmap_prl[, drug] - middle)
  
  cmap_prl %>%
    group_by(SYMBOL) %>%
    arrange(desc(dist)) %>%
    slice(1) %>%
    ungroup() %>%
    select(-dist) %>%
    inner_join(cmap_prls, by="SYMBOL") ->
    cmap_prls
}

class(cmap_prls) <- "data.frame"

for (drug in drug_names) {
  #re-rank (to fill in gaps)
  rank <- order(cmap_prls[, drug])
  cmap_prls[rank, drug] <- 1:nrow(cmap_prls)
}

row.names(cmap_prls) <- cmap_prls[, "SYMBOL"]
cmap_prls <- cmap_prls[, drug_names]
cmap_prls <- as.matrix(cmap_prls)

devtools::use_data(cmap_prls)