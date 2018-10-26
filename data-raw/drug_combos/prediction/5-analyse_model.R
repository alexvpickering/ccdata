library(xgboost)
library(data.table)

setwd("/home/alex/Documents/Batcave/GEO/ccdata/data-raw/drug_combos/prediction")

# load data to train stacker on
train <- fread("data/stack/train.csv")
names(train) <- c("net_preds", "drug1_dprime", "drug2_dprime", "drug1_vardprime",
                  "drug2_vardprime", "combo_dprime")

X <- train[, !"combo_dprime", with=FALSE]

y <- train$combo_dprime

history <- readRDS("data/history.rds")

#preds data.frame
pdf <- data.frame(xgb = history$pred,
                  avg = (train$drug1_dprime + train$drug2_dprime) / 2,
                  actual = y)

#line of best fit for actual vs xgb (and scaled preds)
fitpr <- lm(xgb ~ actual , data=pdf)



#plot actual vs predicted
library(ggplot2)

ggplot(pdf) +
    #scaled vs actual hexbins
    geom_hex(aes(x = xgb, y = actual)) +
    #xgb vs actual line
    geom_abline(intercept = coef(fitpr)[1], slope = coef(fitpr)[2], colour = "red") +
    #perfect line (x = y)
    geom_abline(intercept = 0, slope = 1, colour = "green")


#plot residuals
pdf$resid_pr <- pdf$actual - pdf$xgb

ggplot(pdf) + geom_hex(aes(x = xgb, y = resid_pr))

#plot prediction accuracy as a function of absolute effect size
cuts <- seq(0, 6.1, 0.1)  # 61 bins

get_accuracy <- function(preds, y, cuts) {

    acs <- c()
    for (i in 1:61) {

        #define bin range
        bin <- abs(y) > cuts[i] & abs(y) <= cuts[i+1]

        #get accuracy in bin
        ac <- sum(sign(y[bin]) == sign(preds[bin])) / sum(bin)
        acs <- c(acs, ac)
    }
    return(acs)
}

xgb_accuracy <- get_accuracy(pdf$xgb, y, cuts)
avg_accuracy <- get_accuracy(pdf$avg, y, cuts)


adf <- data.frame(accuracy = c(xgb_accuracy, avg_accuracy),
                  Model = c(rep(c("NNet + XGB", "Avg"), each = length(avg_accuracy))),
                  es = rep(cuts[-62], 2))

adf$Model <- relevel(adf$Model, "Avg")




ap <- ggplot(adf) +
    geom_line(aes(x = es, y = accuracy, colour = Model)) +
    ylab("Classification Accuracy\n") +
    xlab("\nAbsolute Effect Size") +
    scale_x_continuous(breaks = 0:6) +
    scale_y_continuous(breaks = seq(0.5, 1, 0.05)) +
    scale_color_manual(values = c("#E41A1C", "#377EB8")) +
    theme(panel.background = element_rect(fill = "#f8f8f8"),
          plot.background  = element_rect(fill = "#f8f8f8", colour = "#f8f8f8"),
          legend.background = element_rect(fill = "#f8f8f8"),
          panel.grid.major = element_line(colour = "#dddddd"),
          panel.grid.minor = element_line(colour = "#dddddd"))

# save at 1400 width png
ggsave("accuracy.svg", width = 10, height = 5.233, bg="#f8f8f8")

# overall accuracy
sum(sign(y) == sign(pdf$xgb)) / length(y)
sum(sign(y) == sign(pdf$avg)) / length(y)

# 0.8018
# 0.7896 for 'avg' model

#
get_spears <- function(y, yhat) {

    cors <- c()
    for (i in 1:257) {
        a <- (i*11525) - 11524
        b <- (i*11525)
        cors <- c(cors, cor(y[a:b], yhat[a:b], method='spearman'))
    }
    return(cors)
}

get_errors <- function(y, yhat) {

    errors <- c()
    for (i in 1:257) {
        a <- (i*11525) - 11524
        b <- (i*11525)
        errors <- c(errors, sum((sign(y[a:b]) != sign(yhat[a:b]))/length(a:b)))
    }
    return(errors)
}



# errors for each sample
xgb_errors <- get_errors(pdf$xgb, pdf$actual)
median(xgb_errors)
# 0.1901

avg_errors <- get_errors(pdf$avg, pdf$actual)
median(avg_errors)
# 0.1992

edf <- data.frame(value = c(xgb_errors, avg_errors),
                  Model = c(rep("NNet + XGB", length(xgb_errors)),
                            rep("Avg", length(avg_errors))))

ep <- ggplot(edf, aes(x=value, fill = Model)) +
    geom_density(alpha=0.3, show.legend = FALSE) +
    scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
    scale_y_continuous(limits=c(0, 7)) +
    labs(x="\nError Rate", y= "Density\n") +
    theme(panel.background = element_rect(fill = "#f8f8f8"),
          plot.background  = element_rect(fill = "#f8f8f8", colour = "#f8f8f8"),
          legend.background = element_rect(fill = "#f8f8f8"),
          panel.grid.major = element_line(colour = "#dddddd"))



# spearman correlations for each sample
xgb_spears <- get_spears(pdf$xgb, pdf$actual)
median(xgb_spears)
# 0.8153

avg_spears <- get_spears(pdf$avg, pdf$actual)
median(avg_spears)
# 0.8043


sdf <- data.frame(value = c(xgb_spears, avg_spears),
                  Model = rep(c("NNet + XGB", "Avg"),
                              each = length(xgb_spears)))

sp <- ggplot(sdf, aes(x=value, fill = Model)) +
    geom_density(alpha=0.3) +
    scale_fill_manual(values = c("#E41A1C", "#377EB8")) +
    scale_y_continuous(limits=c(0, 7), labels=NULL) +
    labs(x="\nSpearman Correlation", y= "\n") +
    theme(panel.background = element_rect(fill = "#f8f8f8"),
          plot.background  = element_rect(fill = "#f8f8f8", colour = "#f8f8f8"),
          legend.background = element_rect(fill = "#f8f8f8"),
          panel.grid.major = element_line(colour = "#dddddd"))


library(gridExtra)
library(avpick)

# get legend then remove
legend <- get_legend(sp)
sp  <- sp  + theme(legend.position="none")

mult <- arrangeGrob(ep, sp, legend, ncol=3, widths = c(3, 3, 1))

# save at 1400 width png
ggsave("prediction.svg", mult, width = 10, height = 5.625, bg="#f8f8f8")

# save at 1260 width png
ggsave("prediction_thumb.svg", mult, width = 10, height = 5.233, bg="#f8f8f8")
