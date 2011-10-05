# Plot heatmaps of p-values for enrichment across celltypes and states
# Author: Abhishek Sarkar <aksarkar@mit.edu>
library(ggplot2)
library(Cairo)

args <- commandArgs(TRUE)
D <- read.csv(args[1])
D$cell_type <- factor(D$cell_type,
                      levels=c('gm12878', 'h1', 'hepg2', 'hmec', 'hsmm',
                        'huvec', 'k562', 'nhek', 'nhlf'),
                      labels=c('GM12878', 'H1', 'HepG2', 'HMEC', 'HSMM',
                        'Huvec', 'K562', 'NHEK', 'NHLF'))
if (! 'fold' %in% colnames(D)) {
  D$fold <- factor(D$p < .05 / 63, labels=c('', '\U2217'));
}
hi <- 12
p <- (qplot(data=D, x=state, y=cell_type, fill=fold,
            label=sprintf("%.1g", p), size=I(3), 
            geom=c('tile', 'text'), xlab = 'Enhancer state', ylab='Cell type') +
      scale_fill_gradient2(low='#ff0000', mid='white', high=muted('blue'),
                           midpoint=1, limits=c(0, hi)) +
      scale_x_discrete(limits=c('strong', 'weak', 'all')) +
      coord_equal() +
      theme_bw())
CairoPDF(file='out.pdf', width=4, height=6)
print(p)
dev.off()
