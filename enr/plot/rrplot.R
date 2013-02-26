# RR plots for enrichment down the rank list
# Author: Abhishek Sarkar <aksarkar@mit.edu>
library(directlabels)
library(ggplot2)
library(plyr)
library(Cairo)

filter <- function(d, ...) {
  return(head(d[order(d$y, decreasing=TRUE),], n=10))
}

my.bumpup <- function (d, ...) {
  if (nrow(d) > 1) {
    return bumpup(d)
  }
  else {
    return(d)
  }
}

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
  return(ggplot(X, aes(x=total, y=dev, color=celltype)) +
         geom_line(size=I(.25)) +
         geom_hline(yintercept=0, color='black') +
         geom_dl(aes(label=celltype),
                 method=list(cex=.35, 'last.points', 'filter', 'my.bumpup')) +
         scale_x_continuous(name='SNPs ranked by beta',
                            limits=c(0, 200000),
                            breaks=c(0, 50000, 100000, 150000)) +
         scale_y_continuous(name='Cumulative deviation from expected count') +
         theme_bw() +
         theme(strip.background=element_blank(),
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
Cairo(type='pdf', file='rrplot.pdf', dpi=96, width=4 * length(table(E$feature)),
      height=4 * length(table(E$phenotype)), units='in')
print(p)
warnings()
dev.off()
