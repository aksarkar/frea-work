library(ggplot2)
library(Cairo)

source("~/code/enr/plot/theme_nature.R")

args <- commandArgs(TRUE)
data <- read.table(args[1], header=FALSE)
data$fgroup <- sub("_.*", "", data$V9)
data$fgroup <- reorder(data$fgroup, data$fgroup, length)
p <- (ggplot(data, aes(x=fgroup)) +
      stat_bin(aes(y=..count..), geom='point') +
      theme_nature)
pdf()
print(p)
dev.off()
