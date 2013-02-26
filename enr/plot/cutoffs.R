library(directlabels)
library(ggplot2)
library(plyr)
library(scales)

filter <- function(d, ...) {
  return(head(d[order(d$y, decreasing=TRUE),], n=10))
}

my.bumpup <- function (d, ...) {
  if (nrow(d) > 1) {
    return(bumpup(d))
  }
  else {
    return(d)
  }
}

my.log10_trans <- function() {
  trans <- function(d) {
    laply(d, function(x) {
      stopifnot(x >= 0)
      if (x > 0) {
        return(log10(x))
      }
      else {
        return(0)
      }
    })
  }
  inverse <- function(d) {
    laply(d, function(x) 10 ^ x)
  }
  trans_new('my.log10', trans, inverse)
}


args <- commandArgs(TRUE)
d <- read.table(args[1], header=FALSE, sep=' ')
w <- 4 * length(table(d$V2))
d$V2 <- factor(d$V2, levels=c('enh', 'distal_enh', 'dhs+enh', 'dhs', 'distal_dhs', 'dgf', 'distal_dgf', 'faire-seq', 'distal_faire-seq', 'tss+2kb'))
cutoffs <- unique(d$V4)
if (length(args) > 1) {
  ## Plot enrichment of disjoint intervals between ticks
  offset <- 10 ^ (-.5 * log10(2))
} else {
  offset <- 1
}
p <- (ggplot(d, aes(x=offset * V4, y=V6, color=V3)) +
      geom_line(size=I(.25)) +
      geom_hline(yintercept=1, size=I(.25)) +
      geom_dl(aes(label=V3),
              method=list(cex=.35, 'first.points', 'filter', 'my.bumpup')) +
      scale_x_continuous(name='SNP rank cutoff', breaks=cutoffs, labels=cutoffs,
                         trans='my.log10') +
      expand_limits(x=0) +
      scale_y_continuous(name='Fold enrichment') +
      facet_grid(V1 ~ V2, scales='free_y') +
      theme_bw() +
      theme(axis.text.x=element_text(angle=90, hjust=1),
            panel.grid.minor=element_blank(),
            strip.background=element_blank(),
            legend.position='none'))
pdf(height=4 * length(table(d$V1)), width=w)
print(p)
dev.off()
