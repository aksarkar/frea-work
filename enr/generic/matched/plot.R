library(Cairo)
library(directlabels)
library(ggplot2)
library(plyr)
library(scales)

source("~/code/enr/plot/color_roadmap.R")
source("~/code/enr/plot/theme_nature.R")

bumpdown <- function(d,...) {
    if (nrow(d) > 1) {
        d <- calc.boxes(d)[order(d$y, decreasing=TRUE),]
        "%between%" <- function(v,lims)lims[1]<v&v<lims[2]
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
            d$x <- d$x + .1
        }
    }
    d
}

filter <- function(d, ...) {
  return(head(d[order(d$y, decreasing=TRUE),], n=10))
}

cutoffs.plot <- function(X, ylab) {
    ## X$x <- -log10(1 - pchisq(X$V4, 1))
    X$x <- X$V4
    breaks <- unique(X$x)
    (ggplot(X, aes(x=x, y=y, color=V3)) +
     geom_line(size=I(.35 / ggplot2:::.pt)) +
     geom_dl(aes(label=V3),
             method=list(cex=.6, 'last.points', 'filter', 'bumpdown')) +
     scale_x_continuous(name='Negative log association p-value', breaks=breaks, labels=function(x) {sprintf("%.1f", x)}, expand=c(0, 0)) +
     scale_y_continuous(name=ylab, expand=c(0, 0)) +
     scale_color_roadmap +
     theme_nature +
     theme(legend.position='none',
           plot.margin=unit(c(0, margin, 0, 0), "mm"),
           axis.text.x=element_text(angle=90, hjust=1, vjust=.5)
         ))
}

z.plot <- function(X) {
    X$y <- (X$V5 - X$V6) / sqrt(X$V7)
    cutoffs.plot(X, ylab="z-score") + geom_hline(yintercept=1)
}

fold.plot <- function(X) {
    X$y <- X$V5 / X$V6
    cutoffs.plot(X, ylab="Fold enrichment") + geom_hline(yintercept=1)
}

hyperg.plot <- function(X) {
    X$y <- X$V5
    cutoffs.plot(X, ylab="Negative log enrichment p-value") + theme(legend.position='none')
}

slice.plot <- function(X) {
    X$y <- (X$V5 - X$V6) / sqrt(X$V7)
    p <- stats::p.adjust(1 - pnorm(X$y), 'fdr')
    thresh <- max(X$y[which(p > 1e-4)])
    (ggplot(X, aes(x=x, y=y, fill=x)) +
     geom_bar(stat='identity') +
     geom_hline(yintercept=thresh, color='red', size=.35 / ggplot2:::.pt) +
     scale_y_continuous(name="Enrichment z-score") +
     scale_fill_roadmap +
     theme_nature +
     theme(legend.position='none',
           axis.text.x=element_blank()))
}

cell_plot <- function(X) {
  X$x <- factor(X$V3, levels=eid_ordering)
  slice.plot(X) + scale_x_discrete(name="Reference epigenome")     
}

clusters_plot <- function(X) {
  X$x <- factor(X$V3, levels=enh_cluster_ordering)
  slice.plot(X) + scale_x_discrete(name="Enhancer module")
}


args <- commandArgs(TRUE)
X <- read.table(args[1], header=FALSE, sep=' ')
h <- 50
w <- 183
P <- do.call(args[2], list(X))
Cairo(type='pdf', file=gsub("in$", "pdf", args[1]), height=h, width=w, units='mm', dpi='auto')
print(P)
dev.off()
