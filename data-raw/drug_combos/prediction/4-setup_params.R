library(data.table)

setwd("/home/alex/Documents/Batcave/GEO/ccdata/data-raw/drug_combos/prediction")

# load net 1 weights
net1_W1 <- fread("data/net1/W1.csv")
net1_W2 <- fread("data/net1/W2.csv")

net1_b1 <- fread("data/net1/b1.csv")
net1_b2 <- fread("data/net1/b2.csv")

# load net 2 weights
net2_W1 <- fread("data/net2/W1.csv")
net2_W2 <- fread("data/net2/W2.csv")

net2_b1 <- fread("data/net2/b1.csv")
net2_b2 <- fread("data/net2/b2.csv")


# get gene order (same for X and y)
y <- readRDS("prediction/data/y.rds")
genes <- gsub("combo_", "", colnames(y), fixed = TRUE)

net1 <- list(W1 = net1_W1, W2 = net1_W2, b1 = net1_b1, b2 = net1_b2)
net2 <- list(W1 = net2_W1, W2 = net2_W2, b1 = net2_b1, b2 = net2_b2)


net2 <- lapply(net2, function(x) {
    if (ncol(x) == 1) {
        x <- x[, V1]
    } else {
        x <- as.matrix(x)
        colnames(x) <- NULL
    }
    signif(x, 3)
})

devtools::use_data(net1, net2, genes, compress="xz", overwrite = TRUE)



predict.net <- function(W1, W2, b1, b2, X) {
    z2 <- X %*% W1 + b1
    a2 <- pmax(z2/3, z2)
    return (a2 %*% W2 + b2)
}

predict.nets <- function(params, X) {
predict.net(params$n1W1, params$n1W2,
                      params$n1b1, params$n1b2, X)
}


# test it out!
X <- readRDS("data/X.rds")
params <- readRDS("net_params.rds")

system.time(preds <- predict.nets(params, blah))
