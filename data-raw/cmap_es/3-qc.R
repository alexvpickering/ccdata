
# Group Clustering -------------------------------
library(ccdata)

# load data
rma <- readRDS("~/Documents/Batcave/GEO/ccdata/data-raw/cmap_es/rma_processed.rds")
eset <- rma$eset

data("cmap_es")
drugs  <- colnames(cmap_es)

mod <- rma$ebayes$design
svcols <- grepl("SV", colnames(mod))
mod <- mod[, !svcols]

svobj <- rma$svobj
exprs_sva <- crossmeta:::clean_y(exprs(eset), mod, svobj$sv)

# generate pData for eset
pdata <- data.frame(row.names = sampleNames(eset),
                    stringsAsFactors = FALSE)

pdata$Group <- NA

for (i in 1:1309) {
    samples <- mod[, i] == 1
    pdata[samples, "Group"] <- drugs[i]
}

ctl <- mod[, 1310] == 1
pdata[ctl, "Group"] <- "ctl"
pData(eset) <- pdata



# get mds plot data
tested_high <- c("tretinoin", "LY-294002", "sirolimus", "resveratrol", "vorinostat")
tested_low  <- c("dexamethasone", "doxorubicin", "metformin", "progesterone", "quercetin")

get_mds <- function(tested, eset, exprs_sva) {
    
    mds_data <- list()
    
    for (drug in tested) {
        filt <- pdata$Group == drug
        
        # take a sample of controls (too many to plot)
        set.seed(0)
        ctl <- sample(row.names(pdata)[pdata$Group == "ctl"], 100)
        ctl <- sampleNames(eset) %in% ctl
        
        eset2 <- eset[, filt | ctl]
        exprs_sva2 <- exprs_sva[, filt | ctl]
        
        # MDS plot
        Group <- factor(pData(eset2)$Group)
        palette <- RColorBrewer::brewer.pal(12, "Paired")
        colours <- palette[Group]
        names(colours) <- Group
        
        # plot MDS
        mds <- limma::plotMDS(exprs_sva2, pch = 19, main = "CMAP", col = colours)
        mds <- data.frame(x=mds$x, y=mds$y, Group=Group, query_drug=drug)
        mds_data[[drug]] <- mds
    }
    return(mds_data)
}


mds_high <- get_mds(tested_high, eset, exprs_sva)
mds_low  <- get_mds(tested_low, eset, exprs_sva)

# put into one df
df <- Reduce(rbind, c(mds_high, mds_low))
df$Group <- ifelse(df$Group == "ctl", "vehicle", "treatment")
df$Group <- relevel(as.factor(df$Group), "vehicle")


library(ggplot2)
pl <- ggplot(df, aes(x, y)) + 
    geom_point(aes(colour=Group,  alpha=Group)) +
    scale_alpha_discrete(range = c(0.3, 0.8)) +
    facet_wrap(~query_drug, ncol=5) +
    scale_colour_manual(values=c("#377EB8", "#E41A1C")) +
    labs(x="\nLeading logFC dim 1", y= "Leading logFC dim 2\n") +
    scale_x_continuous() + 
    theme(panel.background = element_rect(fill = "#f8f8f8"),
          plot.background  = element_rect(fill = "#f8f8f8", colour = "#f8f8f8"),
          panel.grid.major = element_line(colour = "#dddddd"),
          plot.margin = unit(c(0.3,1,0.3,0.3), "cm"),
          legend.position = "none")

ggsave("mds.svg", pl, width=10, height=5, bg="#f8f8f8")

# P-Values -------------------------------

# load data
library(ccdata)
data(cmap_tables1)
data(cmap_tables2)
cmap_tables <- c(cmap_tables1, cmap_tables2)
rm(cmap_tables1, cmap_tables2)

