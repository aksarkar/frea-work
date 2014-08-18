## RR plots for enrichment down the rank list
## Author: Abhishek Sarkar <aksarkar@mit.edu>
library(directlabels)
library(ggplot2)
library(grid)
library(plyr)
library(scales)
library(Cairo)

source('~/code/enr/plot/color_roadmap.R')
source('~/code/enr/plot/theme_nature.R')

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

my.bumpup <- function(d, ...) {
    if (nrow(d) > 1) {
        return(bumpup(d))
    }
    else {
        return(d)
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

fold <- function(X) {
    transform(X, y=count / (1 + expected))
}

z.score <- function(X, Y, f) {
    transform(merge(f(X),
                    ddply(f(Y), .(total, phenotype, celltype, feature),
                          function(x) c(m=mean(x$y), s=sd(x$y)))),
              y=(y - m) / s)
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

rrplot <- function(X, name, zero, cutoff, labels) {
    stopifnot(is.data.frame(X))
    stopifnot(is.numeric(zero))
    stopifnot(cutoff > 0)
    my.margin <- unit(c(7 / ggplot2:::.pt, margin, 0, 0), 'mm')
    X <- X[X$total <= cutoff, ]
    Y <- X[X$total == cutoff, ]
    write.table(rev(Y[order(Y$y),]$celltype), sub('.in.gz$', '.txt', args[1]),
                col.names=FALSE, row.names=FALSE, quote=FALSE)
    P <- (ggplot(X, aes(x=total, y=y, color=factor(celltype))) +
          geom_line(size=I(.35 / ggplot2:::.pt)) +
          geom_hline(yintercept=zero, color='black', size=I(.5 / ggplot2:::.pt)) +
          labels +
          annotate('text', x=-Inf, y=Inf, size=(5 / ggplot2:::.pt), hjust=1, vjust=0,
                   label=sprintf('(10^{%d})', floor(log10(max(X$y)))), parse=TRUE) +
          scale_x_continuous(name='SNP rank by increasing p-value',
                             labels=comma,
                             limits=c(0, cutoff),
                             breaks=seq(0, cutoff, cutoff / 4),
                             expand=c(0, 0)) +
          scale_y_continuous(name=name, labels=function(x) {sub('e.*', '', sprintf('%.1e', x))},
                             expand=c(0, 0)) +
          scale_color_roadmap +
          theme_nature +
          theme(legend.position='none',
                plot.margin=my.margin))
    t <- ggplot_gtable(ggplot_build(P))
    t$layout$clip[t$layout$name == 'panel'] <- 'off'
    grid.draw(t)
}

rrplot.dev <- function(X, type) {
    labels <- geom_dl(aes(label=celltype),
                      method=list(cex=7/16,
                          last.points,
                          calc.boxes,
                          dl.trans(x = x + .1),
                          top(10),
                          'bumpdown'))
    rrplot.device(X, 1, type)
    rrplot(dev(X), 'Cumulative deviation from expected count', 0, 100000, labels)
    dev.off()
}

rrplot.example <- function(X, type) {
    labels <- geom_dl(aes(label=celltype), method=list(cex=7/16, 'lines2'))
    rrplot.device(X, 1.4, type)
    rrplot(dev(X), 'Cumulative deviation from expected count', 0, 6.2e6)
    dev.off()
}

rrplot.device <- function(X, aspect, type) {
    panelsize <- 89 - margin
    h <- panelsize * length(table(X$phenotype))
    w <- (aspect * panelsize + margin) * length(table(X$feature))
    Cairo(file=sub('.in.gz$', sprintf('.%s', type), args[1]), type=type, width=w, height=h, units='mm', dpi='auto', family='Helvetica')
}

margin <- 10
args <- commandArgs(TRUE)
D <- read.csv(gzfile(args[1]), header=FALSE)
colnames(D) <- c('total', 'phenotype', 'celltype', 'feature', 'count', 'expected')
do.call(args[2], list(D, 'pdf'))
