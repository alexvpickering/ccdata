library(crossmeta)
library(Biobase)
library(limma)
library(metaMA)
gse_name <- 'GSE85871'
data_dir <- 'data-raw/tcm_es'

get_raw(gse_name, data_dir)
eset <- load_raw(gse_name, data_dir)[[1]]
pdata <- pData(eset)

trt_names <- gsub('_rep.+?$', '', pdata$title)
trt_namesv  <- make.names(trt_names)
ctrl <- 'MCF7_vehicle(DMSO)'

treatment <- trt_namesv

mod <- model.matrix(~0 + treatment)
colnames(mod)  <- gsub("^treatment", "", colnames(mod))
row.names(mod) <- 1:nrow(mod)

#generate null model matrix for SVA
mod0 <- model.matrix(~1, data=pdata)


# surrogate variable analysis
svobj <- sva::sva(exprs(eset), mod, mod0)

#add SVs to mod
modsv <- cbind(mod, svobj$sv)
colnames(modsv) <- c(colnames(mod), paste("SV", 1:svobj$n.sv, sep=""))

#generate contrast names (must be "valid")
# contrasts <- paste(drugs_val, "ctl", sep="-")
trt_names <- unique(trt_names)
trt_names <- setdiff(trt_names, ctrl)
trt_val  <- make.names(trt_names)
contrasts <- paste(trt_val, make.names(ctrl), sep="-")


#run limma analysis
ebayes_sv <- crossmeta:::fit_ebayes(eset, contrasts, modsv)


# combine
df <- ebayes_sv$df.residual + ebayes_sv$df.prior
cmap_tables <- list()

for (i in seq_along(contrasts)) {
    cat('Working on', i, 'of', length(contrasts), '\n')
    drugv <- gsub('-.+$', "", contrasts[i])
    ctrlv <- gsub('^.+-', "", contrasts[i])

    drug <- trt_names[i]


    #get top table
    top_table <- topTable(ebayes_sv, coef=i, n=Inf)

    #add dprime and vardprime
    ni <- sum(mod[, ctrlv])
    nj <- sum(mod[, drugv])
    top_table[,c("dprime", "vardprime")] <- effectsize(top_table$t, ((ni * nj)/(ni + nj)), df)[, c("dprime", "vardprime")]

    #store (use eset probe order)
    cmap_tables[[drug]] <- top_table[featureNames(eset), ]
}

# tcm_es

# cmap_es --------------------------------------

#get dprimes and adjusted p-values
es_probes <- lapply(cmap_tables, function(x) x[, c("adj.P.Val", "dprime")])
es_probes <- do.call(cbind, es_probes)


#add symbol
es_probes[,"SYMBOL"] <- fData(eset)$SYMBOL

# where symbol duplicated, keep smallest p-value
es_probes <- as.data.table(es_probes)
dp   <- grep("dprime$", names(es_probes), value = TRUE)
pval <- grep("adj.P.Val$", names(es_probes), value = TRUE)

tcm_es <- es_probes[, Map(`[`,
                           mget(dp),
                           lapply(mget(pval), which.min)),
                     by = SYMBOL]

tcm_es <- tcm_es[!is.na(tcm_es$SYMBOL), ]


# use symbol for row names
class(tcm_es) <- "data.frame"
row.names(tcm_es) <- tcm_es$SYMBOL

# remove dprime from column names
colnames(tcm_es) <- gsub(".dprime", "", colnames(tcm_es))
tcm_es <- as.matrix(tcm_es[, -1])

# save results
saveRDS(tcm_es, 'data-raw/tcm_es/tcm_es.rds')
