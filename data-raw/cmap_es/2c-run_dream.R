args <- commandArgs(trailingOnly=TRUE)
part <- as.numeric(args)

library(variancePartition)
library(BiocParallel)
param <- SerialParam()

# destructure
a <- readRDS('cmap_es/dream/dream_args.rds')
form <- a$form
pdata <- a$pdata
L <- a$L
colnames(L) <- paste0('L', seq_len(ncol(L)))
Linit <- a$Linit
thetaInit <- a$fitInit@theta
fixefInit <- lme4::fixef(a$fitInit)
rm(a); gc()

# run as parts of size 300 (~75 total)
it <- seq(1, 22268)
init <- (part-1) * 300
iend <- min(22268, init + 299)
it <- it[init:iend]

rpath <- '/n/scratch2/ap491/ccdata/data-raw/cmap_es/dream/resLists'

for (i in seq_along(it)) {
  cat('Working on', it[i], 'of', tail(it, 1), '...\n')

  fpath <- file.path(rpath, paste0(it[i], '.rds'))
  exprs_fpath <- file.path(rpath, paste(it[i], 'exprs.rds', sep = '_'))
  if (file.exists(fpath)) next

  # get next gene
  exprs <- readRDS(exprs_fpath)
  resl  <- dream(exprs, form, pdata, L = L,
                 Linit = Linit,
                 thetaInit = thetaInit,
                 fixefInit = fixefInit,
                 return.resList = TRUE,
                 BPPARAM = param)

  # save
  saveRDS(resl, fpath)
  rm(resl, exprs); gc()
}