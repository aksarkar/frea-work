# Plot heatmaps of p-values for enrichment across celltypes and states
# Author: Abhishek Sarkar <aksarkar@mit.edu>
library(ggplot2)
library(Cairo)

CairoFonts(regular="DejaVu Sans")
args <- commandArgs(TRUE)
D <- read.csv(args[1])
D$external <- factor(D$state,
                     levels=c('11','12','11,12','4','5','4,5','4,5,11,12'),
                     labels=c('4', '5', 'strong', '6', '7', 'weak', 'all'))
D$cell_type <- factor(D$cell_type,
                      levels=c('GM12878', 'H1', 'HepG2', 'HMEC', 'HSMM',
                        'Huvec', 'K562', 'NHEK', 'NHLF'))
if (! 'fold' %in% colnames(D)) {
  D$fold <- factor(D$p < .05 / 63, labels=c('', '\U2217'));
}
hi <- ceiling(max(-log10(D$p)))
p <- (qplot(data=D, x=external, y=cell_type, fill=-log10(p), label=fold, size=I(3), 
            geom=c('tile', 'text'), xlab = 'Enhancer state', ylab='Cell type') +
      scale_fill_gradient(low='white', high='blue', limits=c(0, hi)) +
      scale_x_discrete(limits=c('strong', 'weak', 'all')) +
      coord_equal() +
      theme_bw())
Cairo(file='out.pdf', type='pdf', dpi=96, width=4, height=6, units='in')
print(p)
dev.off()
