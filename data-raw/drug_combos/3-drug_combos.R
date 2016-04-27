library(RSQLite)
library(data.table)
library(xgboost)
library(ccmap)
library(ccdata)


# Setup -------------------------

#load model & cmap_full data
data(combo_model)
data(cmap_full)

data   <- cmap_full$data
drugs  <- cmap_full$drugs
probes <- row.names(data)

#indices for drugs
inds <- sapply(seq_along(drugs), function(x) c(x*7-6, x*7))
colnames(inds) <- drugs

#list of unique drug combos
combo_pairs <- combn(drugs, 2, simplify=F)
combo_names <- sapply(combo_pairs, function(x) paste(x[1], x[2], sep=" + "))
names(combo_pairs) <- combo_names

#load annotation information
suppressMessages(library("hgu133a.db"))
suppressMessages(map <- AnnotationDbi::select(hgu133a.db, probes, "SYMBOL"))
map <- map[!is.na(map$SYMBOL),]
genes <- unique(map$SYMBOL)

#make blank table in SQLite database
db.pdc <- dbConnect(SQLite(), dbname="genes_drug_combos.sqlite")
cols <- paste(shQuote(genes), "FLOAT", collapse=", ")

statement <- paste("CREATE TABLE combo_preds",
                                     "(drug_combo TEXT, ", cols, ")", sep="")

dbSendQuery(db.pdc, statement)


pb  <- txtProgressBar(min=1, max=654, style=3)

# 856086 combo_pairs = 1309 * 654
for (i in 1:654) {

    # Probe Predictions -------------------------


    #get array data for 1309 combos
    pairs <- combo_pairs[(i*1309-1308):(i*1309)]
    X     <- combine_cmap_pairs(pairs, data, inds)

    #make predictions with combo_model
    probe_preds <- predict(combo_model, X) - 0.5

    #swap names of drug1 and drug2 columns
    d1_cols <- grep("drug1", colnames(X), value=T)
    d2_cols <- grep("drug2", colnames(X), value=T)
    colnames(X) <- c(d2_cols, d1_cols)

    #restore ordering of drug1 and drug2 columns
    X <- X[, c(d1_cols, d2_cols)]

    #make predictions with combo_model
    probe_preds_rev <- predict(combo_model, X) - 0.5

    #get average & cleanup
    probe_preds <- (probe_preds + probe_preds_rev) / 2
    rm(X, d1_cols, d2_cols, probe_preds_rev)
    gc()

    #probe_preds vector to df
    dim(probe_preds) <- c(length(probes), length(pairs))
    colnames(probe_preds)  <- names(pairs)
    row.names(probe_preds) <- probes
    probe_preds <- as.data.frame(probe_preds)


    # Probes to Genes -------------------------


    #annotate probe probe_preds with SYMBOL
    probe_preds <- probe_preds[map$PROBEID, ]
    probe_preds$SYMBOL <- map$SYMBOL

    #where duplicated SYMBOL, choose probe with prediction furthest from 0.5
    probe_preds <- data.table(probe_preds)
    gene_preds  <- probe_preds[, lapply(.SD, function (col) col[which.max(abs(col))]), by='SYMBOL']



    # Add to SQL Database -------------------------


    #format table
    gene_preds <- as.data.frame(gene_preds)
    row.names(gene_preds) <- gene_preds$SYMBOL
    gene_preds <- gene_preds[, colnames(gene_preds) != "SYMBOL"]

    gene_preds <- as.data.frame(t(gene_preds))
    gene_preds$drug_combo <- row.names(gene_preds)
    row.names(gene_preds) <- NULL

    dbWriteTable(db.pdc, "combo_preds", gene_preds[, c("drug_combo", genes)], append = TRUE)

    rm(gene_preds, probe_preds, pairs, X)
    gc()

    setTxtProgressBar(pb, i)
}



