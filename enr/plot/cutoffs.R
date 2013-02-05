library(ggplot2)
args <- commandArgs(TRUE)
d <- read.table(args[1], header=FALSE, sep=' ')
cutoffs <- unique(d$V4)
# Plot enrichment of disjoint intervals between ticks
offset <- 10 ^ (-.5 * log10(2))
p <- (qplot(data=d, x=V4, y=V6, color=V3,
            xlab='SNP rank cutoff',
            ylab='Fold enrichment', geom='line', size=I(.25), facets=V1 ~ V2) +
      scale_x_log10(breaks=cutoffs, labels=cutoffs) +
      geom_hline(yintercept=1, size=I(.25)) +
      theme_bw() +
      theme(axis.text.x=element_text(angle=90, hjust=1),
            panel.grid.minor=element_blank(),
            strip.background=element_blank(),
            legend.position='none'))
pdf(height=4 * length(table(d$V1)), width=4 * length(table(d$V2)))
print(p)
dev.off()
