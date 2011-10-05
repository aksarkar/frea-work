# Plot heatmaps of p-values for enrichment across celltypes and states
# Author: Abhishek Sarkar <aksarkar@mit.edu>
library(ggplot2)
library(Cairo)

heatmap <- function(df, hi, mid) {
  p <- (qplot(data=df, x=cell_type, y=elim, fill=-log10(p),
              label=sprintf("%.1g", p), size=I(3),
              geom=c('tile', 'text'), xlab='Eliminated cell type', ylab='Cell type') +
        scale_fill_gradient2(high='#ff0000', mid='white', low='#b0b0ff',
                             limits=c(0, hi), midpoint=mid) +
        coord_equal() +
        theme_bw() +
        opts(title=df$state[1]))
  print(p)
}

args <- commandArgs(TRUE)
df <- read.csv(args[1])
df$cell_type <- factor(df$cell_type,
                      levels=c('gm12878', 'h1', 'hepg2', 'hmec', 'hsmm',
                        'huvec', 'k562', 'nhek', 'nhlf'),
                      labels=c('GM12878', 'H1', 'HepG2', 'HMEC', 'HSMM',
                        'Huvec', 'K562', 'NHEK', 'NHLF'))
df$elim <- factor(df$elim,
                  levels=c('nhlf','nhek','k562','huvec','hsmm','hmec','hepg2','h1','gm12878'),
                  labels=c('NHLF','NHEK','K562','Huvec','HSMM','HMEC','HepG2','H1','GM12878'))
hi <- ceiling(max(-log10(df$p)))
mid <- -log10(.05 / 27)

Cairo(file='out.pdf', type='pdf', dpi=96, width=8, height=8, units='in')
d_ply(df, .(state), heatmap, hi, mid)
dev.off()
