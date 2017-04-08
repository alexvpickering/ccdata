library(ggplot2)
library(reshape2)
library(ccmap)
library(ccdata)

#---------------

#load in cmap data
data(cmap_es)
data(atc4)

drugs <- colnames(cmap_es)

#get similarity matrix long
tsim <- get_tan_sim(cmap_es, es=TRUE)

#get tp/fp rates
trates <- get_rates(tsim, atc4)

#put tp/fp rates into dataframe
trates_df <- data.frame(fp=trates$fp, tp=trates$tp)

#plot tp vs fp
ggplot(trates_df, aes(fp, tp)) + geom_point()