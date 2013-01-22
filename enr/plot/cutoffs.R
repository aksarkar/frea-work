library(ggplot2)
args <- commandArgs(TRUE)
d <- read.table(args[1], header=FALSE, sep=' ')
p <- (qplot(data=d, x=-log10(V4), y=V5, color=V3,
            xlab='Negative log-transformed p-value cutoff',
            ylab='Fold enrichment', geom='line') +
      facet_grid(V1 ~ V2) +
      theme_bw() +
      opts(strip.background=theme_rect(fill=NA, colour=NA),
           legend.position='none'))
pdf(height=4 * length(table(d$V1)), width=4 * length(table(d$V2)))
print(p)
dev.off()
