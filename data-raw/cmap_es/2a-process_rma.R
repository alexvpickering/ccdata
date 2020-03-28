library(sva)
library(Biobase)
library(metaMA)
library(crossmeta)
library(data.table)

setwd("~/Documents/Batcave/GEO/ccdata/data-raw/")

# Load Data -------------------------


#load RMA processed data for each platform
ht_hga_ea <- readRDS("cmap_es/rma_HT_HG-U133A_EA.rds")
ht_hga <- readRDS("cmap_es/rma_HT_HG-U133A.rds")
hga <- readRDS("cmap_es/rma_HG-U133A.rds")

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

platform <- c(rep('HG-U133A', ncol(hga)),
              rep('HT_HG-U133A', ncol(ht_hga)),
              rep('HT_HG-U133A_EA', ncol(ht_hga_ea)))

#generate eset from all_exprs
all_exprs <- as.matrix(all_exprs)
eset <- new("ExpressionSet", exprs = all_exprs)

# Setup Analysis -------------------------


#generate model matrix
cmap_instances <- read.table("raw/cmap_instances_02.csv",
                             header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)

# remove ' from scan names
cmap_instances$perturbation_scan_id <- gsub("'", '', cmap_instances$perturbation_scan_id)
cmap_instances$vehicle_scan_id <- gsub("'", '', cmap_instances$vehicle_scan_id)

cmap_instances$batch_id <- gsub("a|b", '', cmap_instances$batch_id)


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
  pdata[pids, "cell"]  <- drug_instances$cell
  pdata[pids, "batch"] <- drug_instances$batch_id

  #get aggregated control ids
  cids_ag <- drug_instances$vehicle_scan_id

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

saveRDS(pdata, 'cmap_es/pdata.rds')
saveRDS(all_exprs, 'cmap_es/all_exprs.rds')

#generate model matrix ---
# treatment <- factor(pdata$drugv, levels = c(drugs_val, 'ctl'))

trt_names   <- paste(pdata$drug, pdata$cell, paste0(pdata$molar, 'M'), paste0(pdata$hours, 'h'), sep='_')
trt_namesv  <- make.names(trt_names)

treatment <- factor(trt_namesv, levels = unique(trt_namesv))
batch <- factor(pdata$batch, levels = unique(pdata$batch))

mod <- model.matrix(~0 + treatment + batch)
colnames(mod)  <- gsub("^treatment", "", colnames(mod))
row.names(mod) <- 1:nrow(mod)

#generate null model matrix for SVA
# mod0 <- model.matrix(~1, data=pdata)

pData(eset) <- pdata

rm(all_exprs, ht_hga_ea, ht_hga, hga); gc()

# Analysis -------------------------


#perform sva
# svobj <- sva(exprs(eset), mod, mod0)

#add SVs to mod
# modsv <- cbind(mod, svobj$sv)
# colnames(modsv) <- c(colnames(mod), paste("SV", 1:svobj$n.sv, sep=""))

#generate contrast names (must be "valid")
# contrasts <- paste(drugs_val, "ctl", sep="-")
trt_names  <- unique(trt_names)
trt_namesv <- unique(trt_namesv)
ctl_names <- gsub("^.+?_(.+?)_.+?_([0-9]+?h)", "ctl_\\1_0M_\\2", trt_namesv)

is_ctrl    <- trt_namesv %in% ctl_names
trt_names  <- trt_names[!is_ctrl]
trt_namesv <- trt_namesv[!is_ctrl]
ctl_names  <- ctl_names[!is_ctrl]

contrasts <- paste(trt_namesv, ctl_names, sep="-")


#run limma analysis (2+ hours)
# ebayes_sv <- crossmeta:::fit_ebayes(eset, contrasts, modsv)
ebayes_sv <- crossmeta:::fit_ebayes(eset, contrasts, mod)


#save results
# rma_processed <- list(eset=eset, svobj=svobj, ebayes_sv=ebayes_sv)
# saveRDS(rma_processed, "cmap_es/rma_processed.rds")
rma_processed <- list(eset=eset, ebayes_sv=ebayes_sv)
saveRDS(rma_processed, "cmap_es/rma_processed_ind.rds")




# Combine -------------------------


#values to calc dprime
# df <- ebayes_sv$df.residual + ebayes_sv$df.prior
# ni <- sum(mod[, "ctl"])

# cmap_tables <- list()
# for (i in seq_along(drugs)) {

#   #get top table
#   top_table <- topTable(ebayes_sv, coef=i, n=Inf)

#   #add dprime and vardprime
#   nj <- sum(mod[, i])
#   top_table[,c("dprime", "vardprime")] <- effectsize(top_table$t, ((ni * nj)/(ni + nj)), df)[, c("dprime", "vardprime")]

#   #store (use eset probe order)
#   cmap_tables[[drugs[i]]] <- top_table[featureNames(eset), ]
# }
df <- ebayes_sv$df.residual + ebayes_sv$df.prior
cmap_tables <- list()

for (i in seq_along(contrasts)) {
  cat('Working on', i, 'of', length(contrasts), '\n')
  drug <- gsub('-.+$', "", contrasts[i])
  ctrl <- gsub('^.+-', "", contrasts[i])

  drugv <- trt_names[i]


  #get top table
  top_table <- topTable(ebayes_sv, coef=i, n=Inf)

  #add dprime and vardprime
  ni <- sum(mod[, ctrl])
  nj <- sum(mod[, drug])
  top_table[,c("dprime", "vardprime")] <- effectsize(top_table$t, ((ni * nj)/(ni + nj)), df)[, c("dprime", "vardprime")]

  #store (use eset probe order)
  cmap_tables[[drug]] <- top_table[featureNames(eset), ]
}

saveRDS(cmap_tables, 'cmap_es/cmap_tables_ind.rds')





# cmap_es --------------------------------------

# get map
ensql <- '/home/alex/Documents/Batcave/GEO/crossmeta/data-raw/entrezdt/ensql.sqlite'
annotation(eset) <- 'GPL96'
fData(eset)$PROBE <- featureNames(eset)
sampleNames(eset) <- paste0('s', 1:ncol(eset))

map <- fData(symbol_annot(eset, ensql = ensql))
map <- map[, c('PROBE', 'SYMBOL')]
map <- map[!is.na(map$SYMBOL), ]

#get dprimes and adjusted p-values
es_probes <- lapply(cmap_tables, function(x) x[, c("adj.P.Val", "dprime")])
es_probes <- do.call(cbind, es_probes)


#add symbol
es_probes <- es_probes[map$PROBE, ] #expands 1:many
es_probes[,"SYMBOL"] <- map$SYMBOL

# where symbol duplicated, keep smallest p-value
es_probes <- as.data.table(es_probes)
dp   <- grep("dprime$", names(es_probes), value = TRUE)
pval <- grep("adj.P.Val$", names(es_probes), value = TRUE)

cmap_es <- es_probes[, Map(`[`,
                           mget(dp),
                           lapply(mget(pval), which.min)),
                     by = SYMBOL]


# use symbol for row names
class(cmap_es) <- "data.frame"
row.names(cmap_es) <- cmap_es$SYMBOL

# remove dprime from column names
colnames(cmap_es) <- gsub(".dprime", "", colnames(cmap_es))
cmap_es <- as.matrix(cmap_es[, -1])
colnames(cmap_es) <- trt_names

# save results
cmap_es <- signif(cmap_es, 5)
saveRDS(cmap_es, 'cmap_es/cmap_es_ind.rds')
#devtools::use_data(cmap_es)



# cmap_var --------------------------------------

#get dprimes and adjusted p-values
var_probes <- lapply(cmap_tables, function(x) x[, c("adj.P.Val", "vardprime")])
var_probes <- do.call(cbind, var_probes)


#add symbol
map <- AnnotationDbi::select(hgu133a.db, row.names(var_probes), "SYMBOL")
map <- map[!is.na(map$SYMBOL), ]
var_probes <- var_probes[map$PROBEID, ] #expands 1:many
var_probes[,"SYMBOL"] <- toupper(map$SYMBOL)

# where symbol duplicated, keep smallest p-value
var_probes <- as.data.table(var_probes)
dp   <- grep("vardprime$", names(var_probes), value = TRUE)
pval <- grep("adj.P.Val$", names(var_probes), value = TRUE)

cmap_var <- var_probes[, Map(`[`,
                             mget(dp),
                             lapply(mget(pval), which.min)),
                       by = SYMBOL]


#use symbol for row names
class(cmap_var) <- "data.frame"
row.names(cmap_var) <- cmap_var$SYMBOL

#remove dprime from column names
colnames(cmap_var) <- gsub(".vardprime", "", colnames(cmap_var))
# cmap_var <- as.matrix(cmap_var[, drugs])
cmap_var <- as.matrix(cmap_var[, trt_names])

#save results
cmap_var <- signif(cmap_var, 5)
saveRDS(cmap_var, 'cmap_es/cmap_var_ind.rds')
#devtools::use_data(cmap_var)
