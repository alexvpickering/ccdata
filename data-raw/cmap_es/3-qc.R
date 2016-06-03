
# Group Clustering -------------------------------

# load data
rma <- readRDS("~/Documents/Batcave/GEO/ccdata/data-raw/cmap_es/rma_processed.rds")
eset <- rma$eset

cmap_es <- ccdata::cmap_es
drugs   <- colnames(cmap_es)

mod <- rma$ebayes$design
svcols <- grepl("SV", colnames(mod))
mod <- mod[, !svcols]

svobj <- rma$svobj
exprs_sva <- crossmeta:::clean_y(exprs(eset), mod, svobj$sv)

# generate pData for eset
pdata <- data.frame(row.names = sampleNames(eset),
                    stringsAsFactors = FALSE)

pdata$group <- NA

for (i in 1:1309) {
    samples <- mod[, i] == 1
    pdata[samples, "group"] <- drugs[i]
}

ctl <- mod[, 1310] == 1
pdata[ctl, "group"] <- "ctl"
pData(eset) <- pdata

# subset eset for analysed drugs
tested <- c("LY-294002", "metformin", "sirolimus", "resveratrol", "vorinostat",
            "quercetin", "tretinoin", "progesterone", "dexamethasone",
            "doxorubicin")

filt <- grepl("progesterone", pdata$group)

# take a sample of controls (too many to plot)
set.seed(0)
ctl <- sample(sampleNames(eset)[pdata$group == "ctl"], 100)
ctl <- sampleNames(eset) %in% ctl

eset2 <- eset[, filt | ctl]
exprs_sva2 <- exprs_sva[, filt | ctl]


# MDS plot
group <- factor(pData(eset2)$group)
palette <- RColorBrewer::brewer.pal(12, "Paired")
colours <- palette[group]
names(colours) <- group

# Add extra space to right of plot area
graphics::par(mar = c(5, 4, 2, 6))

# plot MDS
limma::plotMDS(exprs_sva2, pch = 19, main = "CMAP", col = colours)
graphics::legend("topright", inset = c(-0.4, 0), legend = levels(group),
                 fill = colours[levels(group)], xpd = TRUE, bty = "n", cex = 0.65)


# P-Values -------------------------------

# load data
library(ccdata)
data(cmap_tables1)
data(cmap_tables2)
cmap_tables <- c(cmap_tables1, cmap_tables2)
rm(cmap_tables1, cmap_tables2)

