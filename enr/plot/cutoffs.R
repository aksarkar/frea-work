library(ggplot2)
args <- commandArgs(TRUE)
d <- read.table(args[1], header=FALSE, sep=' ')
p <- (qplot(data=d, x=-log10(V4), y=V5, color=V3,
            xlab='Negative log-transformed p-value cutoff',
            ylab='Fold enrichment', geom='line') +
      facet_grid(V1 ~ V2) +
      theme_bw())
pdf(width=24, height=24)
print(p)
dev.off()
