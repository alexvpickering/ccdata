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
Linit <- a$Linit
Lret <- a$Lret
thetaInit <- a$fitInit@theta
fixefInit <- lme4::fixef(a$fitInit)
sigGStruct <- pbkrtest::get_SigmaG(a$fitInit)$G
rm(a); gc()

# generate top tables for each contrast
cons <- colnames(L)
tests <- row.names(L)[apply(L, 2, function(col) which(col == 1))]
ctrls <- row.names(L)[apply(L, 2, function(col) which(col == -1))]

treats <- gsub('^treatment', '', tests)
cols <- c("dprime", "vardprime")

# run as parts of size 300 (~75 total)
it <- seq(1, 22268)
init <- (part-1) * 300
iend <- min(22268, init + 299)
it <- it[init:iend]

rpath <- '/n/scratch2/ap491/ccdata/data-raw/cmap_es/dream/resLists'

for (i in seq_along(it)) {
  
  cat('Working on', it[i], 'of', tail(it, 1), '...\n')

  
  resl_path <- file.path(rpath, paste0(it[i], '.rds'))
  topt_path <- file.path(rpath, paste0('tt_', it[i], '.rds'))
  exprs_path <- file.path(rpath, paste(it[i], 'exprs.rds', sep = '_'))

  if (file.exists(topt_path)) next
  resl <- readRDS(resl_path)
  exprs <- readRDS(exprs_path)

  ret <- variancePartition:::format.resList(resl, exprs, Lret, sigGStruct, univariateContrasts = FALSE, computeResiduals = FALSE)

  top_tables <- list()
  for (j in seq_along(cons)) {
    con   <- cons[j]
    test  <- tests[j]
    ctrl  <- ctrls[j]
    treat <- treats[j]

    tt <- limma::topTable(ret, coef = con, sort.by = 'P', number = Inf)

    # get study degrees of freedom and group classes
    df <- ret$df.residual[, con]

    # get sample sizes for groups
    ni <- sum(ret$design[, ctrl])
    nj <- sum(ret$design[, test])

    # bind effect size values with top table
    es <- metaMA::effectsize(tt$t, ((ni * nj)/(ni + nj)), df)[, cols, drop = FALSE]
    tt <- cbind(tt, es)

    top_tables[[treat]] <- tt
  }

    # save
  saveRDS(top_tables, topt_path)
  rm(resl, exprs, ret, top_tables); gc()
}

quit(status=0)


# join all in single process
it <- seq(1, 22268)

rpath <- '/n/scratch2/ap491/ccdata/data-raw/cmap_es/dream/resLists'

top_tables <- list()
for (i in seq_along(it)) {
  
  cat('Working on', it[i], 'of', tail(it, 1), '...\n')
  topt_path <- file.path(rpath, paste0('tt_', it[i], '.rds'))
  tts <- readRDS(topt_path)

  for (treat in names(tts)) 
    top_tables[[treat]] <- rbind(top_tables[[treat]], tts[[treat]])
}

# update adjusted p.values for correct number of tests
for (i in seq_along(top_tables)) {
  top_tables[[i]]$adj.P.Value <- p.adjust(top_tables[[i]]$adj.P.Value, method = 'BH')
}

saveRDS(top_tables, 'cmap_es/dream/top_tables.rds')