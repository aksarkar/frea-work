library(ggplot2)
library(scales)
library(directlabels)
library(Cairo)

source('~/code/enr/plot/theme_nature.R')
source('~/code/enr/plot/color_roadmap.R')

top <- function(n) {
  function(d, ...) {
    head(d[order(d$y, decreasing=TRUE),], n=n)
  }
}

my.density <- function(X) {
  X$V1 <- factor(X$V1, levels=enh_cluster_ordering)
  (ggplot(X, aes(x=V2, color=V1)) +
   geom_density(size=.25) +
   labs(y='Density') +
   scale_color_roadmap +
   theme_nature +
   theme(axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5),
         legend.position='none'))
}

daf.plot <- function(X) {
  direct.label(my.density(X) + labs(x='DAF'), method='top.bumptwice')
}

tss.abs.dist.plot <- function(X) {
  (my.density(X) +
   scale_x_continuous(breaks=10 ** seq(1,6), labels=comma) +
   coord_trans(x='log10') +
   labs(x='Absolute distance to closest non-overlapping TSS'))
}

args <- commandArgs(TRUE)
X <- read.table(gzfile(args[1]), header=FALSE, sep=' ')
Cairo(type='pdf', file=sub('.txt.gz$', '.pdf', args[1]), width=89, height=89, units='mm')
do.call(args[2], list(X))
dev.off()
