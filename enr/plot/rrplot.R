# RR plots for enrichment down the rank list
# Author: Abhishek Sarkar <aksarkar@mit.edu>
library(directlabels)
library(ggplot2)
library(grid)
library(plyr)
library(Cairo)

filter <- function(n) {
  function(d, ...) {
    head(d[order(d$y, decreasing=TRUE),], n=n)
  }
}

my.bumpup <- function(d, ...) {
  if (nrow(d) > 1) {
    return(bumpup(d))
  }
  else {
    return(d)
  }
}

dev <- function(X) {
  ddply(X, .(phenotype, feature, celltype), transform,
        y=(count - expected) / max(count))
}

fold <- function(X) {
  transform(X, y=count / (1 + expected))
}

z.score <- function(X, Y, f) {
  transform(merge(f(X),
                  ddply(f(Y), .(total, phenotype, celltype, feature),
                        function(x) c(m=mean(x$y), s=sd(x$y)))),
            y=(y - m) / s)
}

rrplot <- function(X, name, zero, cutoff=30000) {
  stopifnot(is.data.frame(X))
  stopifnot(is.numeric(zero))
  stopifnot(cutoff > 0)
  return(ggplot(X[X$total <= cutoff, ],
                aes(x=total, y=y, color=celltype)) +
         geom_line(size=I(.25)) +
         geom_hline(yintercept=zero, color='black', size=I(.25)) +
         geom_dl(aes(label=celltype),
                 method=list(cex=.8, 'last.points', filter(10), 'my.bumpup')) +
         facet_grid(phenotype ~ feature, scale='free') +
         scale_x_continuous(name='Rank', limits=c(0, 1.6 * cutoff)) +
         scale_y_continuous(name=name) +
         theme_bw() +
         theme(strip.background=element_blank(),
               legend.position="none"))
}

args <- commandArgs(TRUE)
D <- read.csv(args[1], header=FALSE)
colnames(D) <- c('total', 'phenotype', 'celltype', 'feature', 'count', 'expected')
E <- read.csv(args[2], header=FALSE)
colnames(E) <- c(colnames(D), 'rep')
h <- 15
w <- 1.6 * h
Cairo(file=gsub(".in", ".pdf", args[1]), type='pdf', width=w, height=h, units='cm', dpi='auto')
## rrplot(dev(D), 'Normalized deviation', 0)
rrplot(fold(D), 'Fold enrichment', 1)
## rrplot(z.score(D, E, dev), 'Z-score (normalized deviation)', 0)
rrplot(z.score(D, E, fold), 'Z-score (fold enrichment)', 0)
warnings()
dev.off()
