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
    ## X$x <- -pchisq(X$V4, df=1, log.p=TRUE, lower.tail=FALSE)
    X$x <- X$V4
    breaks <- unique(X$x)
    (ggplot(X, aes(x=x, y=y, color=V3)) +
     geom_line(size=I(.25)) +
     geom_dl(aes(label=V3),
             method=list(cex=.6, 'last.points', 'filter', 'bumpdown')) +
     scale_x_continuous(name='Negative log association p-value', breaks=breaks, labels=function(x) {sprintf("%.1f", x)}) +
     scale_y_continuous(name=ylab) +
     scale_color_roadmap +
     expand_limits(x=0) +
     facet_grid(V1 ~ V2, scales='free_y') +
     theme_bw() +
     theme(legend.position='none',
           panel.grid.minor=element_blank(),
           plot.margin=unit(c(0, margin, 0, 0), "in"),
           strip.background=element_blank(),
           axis.text.x=element_text(angle=-90, hjust=0, vjust=.5)
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
    cutoffs.plot(X, ylab="Negative log enrichment p-value")
}

slice.plot <- function(X) {
    X$y <- (X$V5 - X$V6) / sqrt(X$V7)
    X$x <- factor(X$V3, levels=enh_cluster_ordering)
    (ggplot(X, aes(x=x, y=y, color=x)) +
     geom_point(size=I(1)) +
     scale_x_discrete(name="Cluster") +
     scale_y_continuous(name="Enrichment z-score") +
     scale_color_roadmap +
     theme_nature +
     theme(legend.position='none'))
}

args <- commandArgs(TRUE)
X <- read.table(args[1], header=FALSE, sep=' ')
margin <- 1
panelsize <- 4
aspect=1
h <- panelsize * length(table(X$V1))
w <- (aspect * panelsize + margin) * length(table(X$V2))
P <- do.call(args[2], list(X))
t <- ggplot_gtable(ggplot_build(P))
t$layout$clip[t$layout$name == "panel"] <- "off"
pdf(file=gsub("in$", "pdf", args[1]), height=4 * length(table(X$V1)), width=w)
grid.draw(t)
dev.off()
