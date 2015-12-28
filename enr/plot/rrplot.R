## RR plots for enrichment down the rank list
## Author: Abhishek Sarkar <aksarkar@mit.edu>
library(directlabels)
library(ggplot2)
library(grid)
library(plyr)
library(reshape)
library(scales)
library(zoo)
library(Cairo)

source('~/code/enr/plot/color_roadmap.R')
source('~/code/enr/plot/theme_nature.R')

margin <- 7

top <- function(n) {
  function(d, ...) {
    head(d[order(d$y, decreasing=TRUE),], n=n)
  }
}

filter <- function(show) {
  function(d, ...) {
    d[d$celltype %in% show,]
  }
}

bumpdown <- function(d,...) {
  if (nrow(d) > 1) {
    d <- calc.boxes(d)[order(d$y, decreasing=TRUE),]
    d$h <- d$h * 1.1
    d <- calc.borders(d)
    '%between%' <- function(v,lims)lims[1]<v&v<lims[2]
    obox <- function(x,y){
      tocheck <- with(x,c(left,(right-left)/2+left,right))
      tocheck %between% with(y,c(left,right))
    }
    for(i in 2:nrow(d)){
      dif <- d$top[i]-d$bottom[i-1]
      overlap <- c(obox(d[i,],d[i-1,]),obox(d[i-1,],d[i,]))
      if(dif > 0 && any(overlap)){
        d$bottom[i] <- d$bottom[i] - dif
        d$top[i] <- d$top[i] - dif
        d$y[i] <- d$y[i] - dif
      }
    }
  }
  d
}

dev <- function(X) {
  ddply(X, .(phenotype, feature, celltype), transform,
        y=(count - expected) / max(count))
}

derivative <- function(X) {
  ddply(X, .(phenotype, feature, celltype), transform,
        y=c(rollmean(diff(y), 5, align='left'), y[length(y)-4:length(y)]))
}

cumulative.fraction <- function(X) {
  ddply(X, .(phenotype, feature, celltype), transform,
        y=count / max(count))
}

rrplot_annotate <- function(X, cutoff) {
  Y <- X[X$total == cutoff, ]
  write.table(rev(Y[order(Y$y),]$celltype), sub('.in.gz$', '.txt', args[1]),
              col.names=FALSE, row.names=FALSE, quote=FALSE)
}

peaks <- function(series, span=3, ties.method='first') {
  stopifnot((span <- as.integer(span)) %% 2 == 1)
  z <- embed(series, span)
  s <- span %/% 2
  v <- max.col(z, ties.method=ties.method) == 1 + s
  pad <- rep(FALSE, s)
  result <- c(pad, v, pad)
  result
}

rrplot_elbow <- function(X) {
  inflections <- daply(X, .(phenotype, feature, celltype), function(x) {
    hull <- x[order(chull(x$total, x$y)),]
    ddy <- diff(hull$y, differences=2)
    hull[which(diff(sign(ddy)) != 0)[1] + 1,]$total
  })
  rank_cutoff <- max(inflections[!is.na(inflections)])
  print(sprintf("%s %d", args[1], rank_cutoff))
  geom_vline(xintercept=rank_cutoff, color='red', size=I(.5 / ggplot2:::.pt))
}

rrplot <- function(X, zero=0, cutoff=NULL, scales.labels=NULL, direct.labels=NULL, plot_elbow=TRUE) {
  stopifnot(is.data.frame(X))
  stopifnot(is.numeric(zero))
  if (plot_elbow) {
    elbow <- rrplot_elbow(X)
  }
  else {
    elbow <- geom_hline(yintercept=zero, color='black', size=I(0))
  }
  if (is.null(cutoff)) {
      cutoff <- max(X$total)
  }
  else {
      X <- X[X$total <= cutoff, ]
  }
  (ggplot(X, aes(x=total, y=y, color=factor(celltype))) +
   geom_line(size=I(.35 / ggplot2:::.pt)) +
   geom_hline(yintercept=zero, color='black', size=I(.5 / ggplot2:::.pt)) +
   elbow +
   direct.labels +
   scale_x_continuous(labels=comma, limits=c(0, cutoff), breaks=seq(0, cutoff, cutoff / 4), expand=c(0, 0)) +
   scale_y_continuous(expand=c(0, 0)) +
   scale_color_roadmap +
   scales.labels +
   theme_nature +
   theme(legend.position='none',
         axis.text=element_text(size=5),
         plot.margin=unit(c(7 / ggplot2:::.pt, margin, 0, 0), 'mm')))
}

rrplot_dev <- function(X, xlab, cutoff) {
  direct.labels <- geom_dl(aes(label=celltype),
                           method=list(cex=7/16,
                             last.points,
                             calc.boxes,
                             dl.trans(x = x + .1),
                             top(10),
                             'bumpdown'))
  scales.labels <- labs(x=xlab, y='Cumulative deviation')
  rrplot(dev(X), 0, as.numeric(cutoff), scales.labels, direct.labels)
}

rrplot_derivative <- function(X, cutoff=50000) {
  direct.labels <- geom_dl(aes(label=celltype),
                           method=list(cex=5/16,
                             first.points,
                             calc.boxes,
                             top(5),
                             'bumpdown'))
  scales.labels <- labs(x=xlab, y='Cumulative deviation')
  rrplot(derivative(dev(X)), 0, as.numeric(cutoff), labs(x='Top SNPs by p-value', y='Derivative of cumulative deviation'), direct.labels)
}

rrplot_snps <- function(X, cutoff=30000) {
  rrplot_dev(X, 'SNP rank by p-value', cutoff)
}

rrplot_cdf <- function(X) {
    rrplot(cumulative.fraction(X),
           cutoff=30000,
           scales.labels=labs(xlab='SNP rank by p-value', ylab='Cumulative fraction'),
           direct.labels=NULL,
           plot_elbow=FALSE)
}

rrplot_corr <- function(X) {
  scales.labels <- labs(x='WTCCC1 rank by p-value', y='Cumulative fraction recovered')
  (rrplot(cumulative.fraction(X), 0, signif(max(X$total), 2), scales.labels, NULL) +
   geom_abline(slope=1/max(X$total), yintercept=0, linetype='dotted', size=I(.35/ggplot2:::.pt), color='black') +
   guides(color=guide_legend('Threshold')) +
   theme(legend.position='right'))
}

rrplot_loci <- function(X, cutoff=10000) {
  rrplot_dev(X, 'Loci by tag p-value', cutoff)
}

rrplot_example <- function(X) {
  scales.labels <- labs(x='SNP rank by p-value', y='Cumulative deviation')
  rrplot(dev(X), 0, signif(max(X$total), 2), scales.labels, NULL, plot_elbow=FALSE)
}

rrplot_clusters <- function(X) {
  direct.labels <- geom_dl(aes(label=celltype),
                           method=list(cex=5/16,
                             last.points,
                             calc.boxes,
                             dl.trans(x = x + .1),
                             top(10),
                             'bumpdown'))
  Y <- dev(X)
  dY <- ddply(Y, .(phenotype, feature, celltype), transform, dy=c(0, diff(y))^2)
  Z <- cast(dY, celltype ~ total, value='dy')
  h <- hclust(dist(Z))
  groups <- cutree(h, k=8)
  Y$celltype <- groups[Y$celltype]
  W <- ddply(Y, .(total, celltype),
             function(A) {data.frame(y=quantile(A$y, .5),
                                     ymin=quantile(A$y, .25),
                                     ymax=quantile(A$y, .75))})
  (rrplot(W, 0, 100000, labs(x='SNPs by p-value', y='Cumulative deviation'), direct.labels) +
   geom_ribbon(aes(ymin=ymin, ymax=ymax, fill=factor(celltype)), color=NA, alpha=.15) +
   scale_color_brewer(palette='Dark2') +
   scale_fill_brewer(palette='Dark2'))
}

rrplot_device <- function(X, aspect, type) {
  panelsize <- 60 - margin
  h <- panelsize * length(table(X$phenotype))
  w <- (aspect * panelsize + margin) * length(table(X$feature))
  Cairo(file=sub('.in.gz$', sprintf('.%s', type), args[1]), type=type, width=w, height=h, units='mm', dpi='auto', family='Helvetica')
}

logp_ticks <- function(Y) {
  T <- rev(table(cut(Y$V5, labels=seq(1, 10), breaks=seq(0, 10), right=FALSE)))
  Z <- data.frame(thresh=names(T), rank=cumsum(T), row.names=NULL)
  geom_vline(aes(xintercept=rank), data=Z, linetype='dashed', size=I(.35 / ggplot2:::.pt))
}

chisq_ticks <- function(Y) {
  T <- rev(table(cut(Y$V5, labels=seq(1, 10), breaks=qchisq(1 - 10 ** -seq(0, 10), 1), right=FALSE)))
  Z <- data.frame(thresh=names(T), rank=cumsum(T), row.names=NULL)
  ## write.table(Z, file='ticks', col.names=FALSE, row.names=FALSE, quote=FALSEn)
  geom_vline(aes(xintercept=rank), data=Z, linetype='dashed', size=I(.35 / ggplot2:::.pt))
}

id_ticks <- function(Y) {
  colnames(Y) <- c('thresh', 'rank')
  geom_vline(aes(xintercept=rank), data=Y, linetype='dashed', size=I(.35 / ggplot2:::.pt))
}

no_ticks <- function(Y) {
  NULL
}

rrplot_draw <- function(bin.file, assoc.file, aspect, type, ticks.fn, rrplot_fn, ...) {
  X <- read.csv(gzfile(bin.file), header=FALSE, col.names=c('total', 'phenotype', 'celltype', 'feature', 'count', 'expected'))
  Y <- read.table(gzfile(assoc.file), header=FALSE)
  rrplot_device(X, as.numeric(aspect), type)
  P <- do.call(rrplot_fn, list(X, ...)) + do.call(ticks.fn, list(Y))
  t <- ggplot_gtable(ggplot_build(P))
  t$layout$clip[t$layout$name == 'panel'] <- 'off'
  grid.draw(t)
  dev.off()
}

args <- commandArgs(TRUE)
do.call(rrplot_draw, as.list(args))
