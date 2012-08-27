# RR plots for enrichment down the rank list
# Author: Abhishek Sarkar <aksarkar@mit.edu>
library(ggplot2)
library(Cairo)

scale_cell_type <-
  scale_color_manual(name='Cell type',
                     values=c('H1' = '#3e8f54', 'Huvec' = '#e2560f',
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

rrplot <- function(D, colors) {
  return (qplot(data=D, x=total, y=dev, geom='line', size=I(.25), color=feature,
                xlab='SNPs ranked by p-value',
                ylab='Normalized deviation from expected count') + 
          colors +
          scale_x_continuous(limits=c(0, 150000),
                             breaks=c(0, 50000, 100000, 150000)) +
          scale_y_continuous(limits=c(0, max(D$dev))) +
          theme_bw())
}

args <- commandArgs(TRUE)
D <- read.csv(args[1])
D$denom = rep(tail(D$expected, length(table(D$feature))), length.out=length(D$count))
D$dev <- (D$count - D$expected) / D$denom
p <- rrplot(subset(D, total < 150000),
            scale_color_brewer(name='disease x feature x celltype', palette='Set1'))
Cairo(type='pdf', file='out.pdf', dpi=96, width=8, height=5, units='in')
print(p)
dev.off()
