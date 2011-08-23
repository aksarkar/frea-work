# Plot heatmaps of p-values for enrichment across celltypes and states
# Author: Abhishek Sarkar <aksarkar@mit.edu>
library(ggplot2)
library(Cairo)

args <- commandArgs(TRUE)
D <- read.csv(args[1])
hi <- ceiling(max(-log10(D$p)))
p <- (qplot(data=D, x=state, y=cell_type, fill=-log10(p), label=fold, size=I(4),
            geom=c('tile', 'text'), xlab = 'State', ylab='Cell type') +
      scale_fill_gradient(low='white', high='blue', limits=c(0, hi)) +
      coord_equal() +
      theme_bw())
Cairo(file='out.pdf', type='pdf', dpi=96)
print(p)
dev.off()
