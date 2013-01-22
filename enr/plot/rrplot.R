# RR plots for enrichment down the rank list
# Author: Abhishek Sarkar <aksarkar@mit.edu>
library(ggplot2)
library(Cairo)

scale_cell_type <-
  scale_color_manual(name='Cell type',
                     values=c('H1' = '#3e8f54', 'HUVEC' = '#e2560f',
                       'GM12878' = '#9e320e', 'HepG2' = '#1182bb',
                       'HMEC' = '#3fe7d2', 'K562' = '#3c11bc',
                       'NHEK' = '#bb1d94', 'HSMM' = '#e0c31b',
                       'NHLF' = '#806cd3'))

scale_state_aggregate <-
  scale_color_manual(name='State (aggregate)',
                     values=c('promoter' = '#ff0000', 'enhancer' = '#faca00',
                       'insulator' = '#09befe', 'transcribed' = '#00b050',
                       'repressed' = '#7f7f7f', 'other' = 'black',
                       'Poly-A+ RNA-Seq' = '#005c1f', 'diffexpr' = '#0060a0'))

dev <- function(X, binsize=1000) {
  return(data.frame(X$total, (X$count - X$expected) / max(X$count)))
}

rrplot <- function(X) {
  return(qplot(data=X, x=total, y=dev, geom='line', size=I(.25),
               color=celltype, xlab='SNPs ranked by p-value',
               ylab='Normalized deviation from expected count') +
         scale_x_continuous(limits=c(0, 150000),
                            breaks=c(0, 50000, 100000, 150000)) +
         geom_hline(yintercept=0, color='black') +
         theme_bw() +
         opts(strip.background=theme_rect(fill=NA, colour=NA),
              legend.position="none"))
}

args <- commandArgs(TRUE)
D <- read.csv(args[1])
E <- ddply(D, .(phenotype, feature, celltype), dev)
rm(D)
colnames(E) <- c('phenotype', 'feature', 'celltype', 'total', 'dev')
p <- (rrplot(subset(E, total < 150000)) +
      facet_grid(phenotype ~ feature, scale='free'))
write.csv(E, file='rrplot.csv', quote=FALSE, row.names=FALSE)
Cairo(type='pdf', file='out.pdf', dpi=96, width=4 * length(table(E$feature)),
      height=4 * length(table(E$phenotype)), units='in')
print(p)
warnings()
dev.off()
