# Plot heatmaps of p-values for enrichment across celltypes and states
# Author: Abhishek Sarkar <aksarkar@mit.edu>
library(ggplot2)
library(Cairo)

args <- commandArgs(TRUE)
D <- read.csv(args[1])
D$external <- factor(D$state,
                     levels=c('11','12','11,12','4','5','4,5','4,5,11,12'),
                     labels=c('4', '5', 'strong', '6', '7', 'weak', 'all'))
D$cell_type <- factor(D$cell_type,
                      ## levels=c('nhlf', 'nhek', 'k562', 'huvec', 'hsmm', 'hmec', 'hepg2', 'h1', 'gm12878'),
                      levels=c('NHLF', 'NHEK', 'K562', 'Huvec', 'HSMM', 'HMEC', 'HepG2', 'H1', 'GM12878'))
hi <- ceiling(max(-log10(D$p)))
mid <- -log10(.05 / 27)
p <- (qplot(data=D, x=external, y=cell_type, fill=-log10(p),
            label=sprintf("%.1g", p), size=I(3),
            geom=c('tile', 'text'), xlab = 'Enhancer state', ylab='Cell type') +
      scale_fill_gradient2(high='#ff0000', mid='white', low='#b0b0ff',
                           limits=c(0, hi), midpoint=mid) +
      scale_x_discrete(limits=c('strong', 'weak', 'all')) +
      coord_equal() +
      theme_bw())
Cairo(file='out.pdf', type='pdf', dpi=96, width=4, height=6, units='in')
print(p)
dev.off()
