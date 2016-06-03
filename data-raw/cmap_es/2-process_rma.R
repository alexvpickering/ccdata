library(sva)
library(Biobase)
library(metaMA)
library(crossmeta)
library(hgu133a.db)

setwd("~/Documents/Batcave/GEO/ccdata/data-raw/")

# Load Data -------------------------


#load RMA processed data for each platform
ht_hga_ea <- readRDS("cmap-es/rma_HT_HG-U133A_EA.rds")
ht_hga <- readRDS("cmap-es/rma_HT_HG-U133A.rds")
hga <- readRDS("cmap-es/rma_HG-U133A.rds")

#log2 ht_hga (RMA from xps doesn't log2)
exprs(ht_hga) <- log2(exprs(ht_hga))

#fix up sample names
ht_hga_names <- strsplit(sampleNames(ht_hga), "[.]")
ht_hga_names <- sapply(ht_hga_names, function(x) substring(x[1], 2))
ht_hga_names <- gsub("_", ".", ht_hga_names)
sampleNames(ht_hga) <- ht_hga_names

sampleNames(ht_hga_ea) <- gsub(".CEL", "", sampleNames(ht_hga_ea))
sampleNames(hga) <- gsub(".CEL", "", sampleNames(hga))

#merge data from all platforms
all_exprs <- merge(exprs(hga), exprs(ht_hga), by="row.names")
row.names(all_exprs) <- all_exprs$Row.names
all_exprs <- merge(all_exprs, exprs(ht_hga_ea), by="row.names")
row.names(all_exprs) <- all_exprs$Row.names
all_exprs <- all_exprs[,-(1:2)]

#generate eset from all_exprs
all_exprs <- as.matrix(all_exprs)
eset <- new("ExpressionSet", exprs = all_exprs)

# Setup Analysis -------------------------


#generate model matrix
cmap_instances <- read.table("raw/cmap_instances_02.csv",
                             header=TRUE, sep="\t", quote='', fill=TRUE, stringsAsFactors=FALSE)

cmap_instances$batch_id <- (gsub("(\\d+)\\w", "\\1", cmap_instances$batch_id))

#"valid" drug names (needed for limma makeContrasts)
drugs     <- unique(cmap_instances$cmap_name)
drugs_val <- make.names(drugs, unique = TRUE)

pdata <- data.frame(row.names = sampleNames(eset), 
                    drug = character(7056),
                    drugv = character(7056),
                    molar = numeric(7056),
                    hours = numeric(7056),
                    cell = character(7056),
                    batch = numeric(7056),
                    stringsAsFactors = FALSE)


for (i in seq_along(drugs)) {
  drug     <- drugs[i]
  drug_val <- drugs_val[i]

  #get cmap info for drug
  drug_instances <- cmap_instances[cmap_instances$cmap_name == drug, ]

  #add drug cmap info
  pids <- drug_instances$perturbation_scan_id

  pdata[pids, "drug"]  <- drug
  pdata[pids, "drugv"] <- drug_val
  pdata[pids, "molar"] <- drug_instances$concentration..M.
  pdata[pids, "hours"] <- drug_instances$duration..h.
  pdata[pids, "cell"]  <- drug_instances$cell2
  pdata[pids, "batch"] <- drug_instances$batch_id

  #get aggregated control ids
  cids_ag <- drug_instances$vehicle_scan_id4

  #de-aggregate control ids
  for (j in seq_along(cids_ag)) {
    cid_ag <- cids_ag[j]
    pid    <- pids[j]

    cids <- c()
    #if multiple controls: get prefix/sufixes
    if (length(strsplit(cid_ag, "[.]")[[1]]) > 2) {

      pref = strsplit(pid,"[.]")[[1]][1]
      sufs = strsplit(cid_ag,"[.]")
      sufs = sufs[[1]][-which(sufs[[1]] %in% "")]

      #paste prefix/suffixes together
      cids <- c(cids, paste(pref, sufs, sep="."))

      #if single control: add cid directly
    } else {
      cids <- cid_ag
    }
    #add control cmap info
    pdata[cids, "drug"]  <- "ctl"
    pdata[cids, "drugv"] <- "ctl"
    pdata[cids, "molar"] <- rep(0, length(cids))
    pdata[cids, "hours"] <- rep(pdata[pid, "hours"], length(cids))
    pdata[cids, "cell"]  <- rep(pdata[pid, "cell"],  length(cids))
    pdata[cids, "batch"] <- rep(pdata[pid, "batch"], length(cids))
  }
}

#generate model matrix
treatment <- factor(pdata$drugv, levels = c(drugs_val, 'ctl'))

mod <- model.matrix(~0 + treatment + cell)
colnames(mod)  <- gsub("^treatment", "", colnames(mod))
row.names(mod) <- 1:nrow(mod)

#generate null model matrix for SVA
mod0 <- model.matrix(~1, data=pdata)

pData(eset) <- pdata

# Analysis -------------------------


#perform sva
svobj <- sva(exprs(eset), mod, mod0)

#add SVs to mod
modsv <- cbind(mod, svobj$sv)
colnames(modsv) <- c(colnames(mod), paste("SV", 1:svobj$n.sv, sep=""))

#generate contrast names (must be "valid")
contrasts <- paste(drugs_val, "ctl", sep="-")

#run limma analysis (2+ hours)
ebayes_sv <- fit_ebayes(eset, contrasts, modsv)


#save results
rma_processed2 <- list(eset=eset, svobj=svobj, ebayes_sv=ebayes_sv)
saveRDS(rma_processed, "cmap_es/rma_processed.rds")




# Combine -------------------------


#values to calc dprime
df <- ebayes_sv$df.residual + ebayes_sv$df.prior
ni <- sum(mod[, "ctl"])

cmap_tables <- list()
for (i in seq_along(drugs)) {

  #get top table
  top_table <- topTable(ebayes_sv, coef=i, n=Inf)

  #add dprime
  nj <- sum(mod[, i])
  top_table$dprime <- effectsize(top_table$t, ((ni * nj)/(ni + nj)), df)[, "dprime"]

  #store (use eset probe order)
  cmap_tables[[drugs[i]]] <- top_table[featureNames(eset), ]
}
#save
devtools::use_data(cmap_tables)




# Annotate --------------------------------------



#get dprimes and adjusted p-values
es_probes <- lapply(cmap_tables, function(x) x[, c("adj.P.Val", "dprime")])
es_probes <- do.call(cbind, es_probes)


#add symbol
map <- AnnotationDbi::select(hgu133a.db, row.names(es_probes), "SYMBOL")
map <- map[!is.na(map$SYMBOL), ]
es_probes <- es_probes[map$PROBEID, ] #expands 1:many
es_probes[,"SYMBOL"] <- toupper(map$SYMBOL)

# where symbol duplicated, keep smallest p-value
es_probes <- as.data.table(es_probes)
dp   <- grep("dprime$", names(es_probes), value = TRUE)
pval <- grep("adj.P.Val$", names(es_probes), value = TRUE)

cmap_es <- es_probes[, Map(`[`, 
                            mget(dp), 
                            lapply(mget(pval), which.min)),
                      by = SYMBOL]


#use symbol for row names
class(cmap_es) <- "data.frame"
row.names(cmap_es) <- cmap_es$SYMBOL

#remove dprime from column names
colnames(cmap_es) <- gsub(".dprime", "", colnames(cmap_es))
cmap_es <- as.matrix(cmap_es[, drugs])

#save results
devtools::use_data(cmap_es)
