library(variancePartition)
library(BiocParallel)
setwd("~/Documents/Batcave/GEO/ccdata/data-raw/")

pdata <- readRDS('cmap_es/pdata.rds')
all_exprs <- readRDS('cmap_es/all_exprs.rds')

# mixed effect model ---
pdata$platform <- factor(platform)
pdata$cell <- factor(pdata$cell)
pdata$batch <- factor(pdata$batch)
pdata$treatment <- paste(pdata$drug, paste0(pdata$molar, 'M'), paste0(pdata$hours, 'h'), sep='_')

# variance partition (treat all as random effects) ---
form <- ~(1|treatment) + (1|cell) + (1|batch) + (1|platform)
param <- SerialParam()
varPart <- fitExtractVarPartModel(all_exprs, form, pdata, BPPARAM = param)
saveRDS(varPart, 'cmap_es/varPart.rds')
plotVarPart(sortCols(varPart))

# remove effect of batch and platform then redo
varPart <- readRDS('cmap_es/varPart.rds')
param <- SerialParam()
residList <- fitVarPartModel(all_exprs, ~ (1|batch) + (1|platform), pdata, BPPARAM = param, fxn=residuals)
residMatrix <- do.call(rbind, residList)

# fit model on residuals
form <- ~ (1|treatment) + (1|cell)
varPartResid <- fitExtractVarPartModel(residMatrix, form, pdata, BPPARAM = param)
saveRDS(varPartResid, 'cmap_es/varPartResid.rds')
plotVarPart(sortCols(varPartResid))

# differential expression analysis (treatment has to be fixed effect) ----
param <- SerialParam()
form <- ~ 0 + treatment + (1|cell) + (1|batch) + (1|platform)

# exclude treatments with colinearity issues (see below)
keep <- row.names(pdata)[!pdata$treatment %in% maxs]
pdata <- pdata[keep, ]
all_exprs <- all_exprs[, keep]

# univariate contrasts (faster to supply)
Linit <- variancePartition:::.getAllUniContrasts(all_exprs, form, pdata, return.Linit = TRUE)

# interested contrasts
levels <- unique(pdata$treatment)
cons <- paste0('treatment', levels[!grepl('^ctl_', levels)])
L <- lapply(cons, function(con) {
  ctrl <- ifelse(grepl('_6h$', con), 'treatmentctl_0M_6h', 'treatmentctl_0M_12h')
  getContrast(all_exprs, form, pdata, c(con, ctrl), L = Linit)
})
L <- do.call(cbind, L)
colnames(L) <- paste0('L', seq_len(ncol(L)))

# initial fit used to speed up subsequent
system.time(fitInit <- dream(all_exprs, form, pdata, L = L, Linit = Linit, return.fitInit = TRUE, BPPARAM=param))
variancePartition:::checkModelStatus(fitInit, showWarnings=TRUE, dream=TRUE, colinearityCutoff=0.999)
#    user  system elapsed
# 252.463   0.390 230

# levels of treatment with very high correlation to ctrl/each other cause colinearity issues
# re-run above excluding maxs
vcor <- colinearityScore(fitInit)
vcor <- attributes(vcor)[[1]]
diag(vcor) <- 0
maxs <- apply(vcor, 2, max)
maxs <- names(maxs)[maxs > 0.999]
maxs <- setdiff(maxs, 'treatmentctl_0M_6h')
maxs <- gsub('^treatment', '', maxs)


Lret <- variancePartition:::format.contrasts(all_exprs, form, pdata, L, Linit = Linit)

# save arguments to run as parts
# see 2-run_dream.R and 2-process_dream.R
# ran on O2
dream_args <- list(form = form, pdata = pdata, L = L, Linit = Linit, fitInit = fitInit, Lret = Lret)
dir.create('cmap_es/dream')
dir.create('cmap_es/dream/resLists')
saveRDS(dream_args, 'cmap_es/dream/dream_args.rds')
saveRDS(all_exprs, 'cmap_es/all_exprs.rds')

# save expression values for each gene as seperate matrix
rpath <- '/n/scratch2/ap491/ccdata/data-raw/cmap_es/dream/resLists'
for (i in seq_len(nrow(all_exprs))) {
  cat('Working on', i, 'of', nrow(all_exprs), '\n')
  exprs_fpath <- file.path(rpath, paste(i, 'exprs.rds', sep = '_'))
  saveRDS(all_exprs[i,,drop=FALSE], exprs_fpath)
}
