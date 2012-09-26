# Plot ROC curves
# Expects CSV with entries (fpr, tpr, label) as first argument
# Author: Abhishek Sarkar
library(Cairo)
library(ggplot2)

args <- commandArgs(TRUE)
d <- read.csv(args[1], header=0)
p <- (ggplot(d, aes(V1, V2, color=factor(V3))) +
      geom_line(size=.25) +
      scale_color_brewer(palette='Set1', name='Dataset') +
      scale_x_continuous(limits=c(0, 1), name='False positive rate') +
      scale_y_continuous(limits=c(0, 1), name='True positive rate') +
      coord_equal() +
      geom_abline(slope=1, intercept=0, color='gray', size=0.1) +
      theme_bw())
CairoPDF('out.pdf')
print(p)
dev.off()
