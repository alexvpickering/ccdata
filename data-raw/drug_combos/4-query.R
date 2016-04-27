library(RSQLite)
library(data.table)
library(foreach)
library(doMC)
library(ccmap)

setwd("~/Documents/Batcave/GEO/2-cmap/CR")

# Setup -------------------------
registerDoMC(6)

# 856086 combos = 1309 * 109 * 6
res_list <- foreach(i=1:6) %dopar% {

  #load CR dprimes
  es2 <- readRDS("~/Documents/Batcave/GEO/1-meta/CR/es2.rds")
  query_genes <- get_dprimes(es2)$meta

  #connect to db
  db <- dbConnect(SQLite(), dbname="~/Documents/Batcave/GEO/ccdata/data-raw/drug_combos/drug_combos.sqlite")

  a <- i*109-108
  b <- i*109

  pb  <- txtProgressBar(min=a, max=b, style=3)
  res <- list()

  for (j in a:b) {

    #get preds for drug combos
    statement   <- paste("SELECT * from combo_preds WHERE rowid BETWEEN", (j*1309)-1308, "AND", j*1309)
    combo_preds <- dbGetQuery(db, statement)

    #here I do some stuff to the result returned from the query
    combo_names <- combo_preds$drug_combo
    combo_preds <- as.data.frame(t(combo_preds[,-1]))

    colnames(combo_preds)  <- combo_names

    #get top drug combos
    top_combos <- get_top_drugs(query_genes, as.matrix(combo_preds))

    #update progress and store result
    setTxtProgressBar(pb, j)
    res[[ length(res)+1 ]] <- top_combos
  }

  dbDisconnect(db)
}