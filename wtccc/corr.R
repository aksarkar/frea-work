library(plyr)
library(ggplot2)
library(Cairo)

source("~/code/enr/plot/theme_nature.R")

args <- commandArgs(TRUE)
data <- read.delim(gzfile(args[1]), header=FALSE)
data$V6 <- order(data$V4, decreasing=TRUE)
data$V7 <- order(data$V5, decreasing=TRUE)
data$top.10.percentile <- data$V7 < .1 * nrow(data)
recovered <- data.frame(x=cumsum(as.numeric(data[data$V6,]$top.10.percentile)))
panelsize <- 63
Cairo(type='pdf', file='plot.pdf', width=panelsize, height=panelsize, units='mm')
(ggplot(recovered, aes(x=x)) +
 stat_ecdf(geom='step') +
 theme_nature)
dev.off()
