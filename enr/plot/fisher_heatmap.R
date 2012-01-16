# Plot heatmaps of p-values for enrichment across celltypes and states
# Author: Abhishek Sarkar <aksarkar@mit.edu>
library(ggplot2)
library(Cairo)

args <- commandArgs(TRUE)
D <- read.csv(args[1])
D$cell_type <- factor(D$cell_type,
                      levels=c('nhlf', 'nhek', 'k562', 'huvec', 'hsmm', 'hmec',
                        'hepg2', 'h1', 'gm12878'),
                      labels=c('NHLF', 'NHEK', 'K562', 'Huvec', 'HSMM', 'HMEC',
                        'HepG2', 'H1', 'GM12878'))
hi <- 12
p <- (qplot(data=D, x=state, y=cell_type, fill=fold,
            label=sprintf("%.1g", p), size=I(3), 
            geom=c('tile', 'text'), xlab = 'Enhancer state', ylab='Cell type') +
      scale_fill_gradient2(high='#ff0000', mid='white', low='#b0b0ff',
                           midpoint=1, limits=c(0, hi)) +
      scale_x_discrete(limits=c('strong', 'weak', 'all')) +
      coord_equal() +
      theme_bw())
CairoPDF(file='out.pdf', width=4, height=6)
print(p)
dev.off()
