library(directlabels)
library(ggplot2)
library(plyr)
library(scales)
source("~/code/enr/plot/color_roadmap.R")

filter <- function(d, ...) {
  return(head(d[order(d$y, decreasing=TRUE),], n=10))
}

cutoffs.plot <- function(X, y, ylab) {
    (ggplot(X, aes(x=x, y=y, color=V3)) +
     geom_line(size=I(.25)) +
     geom_dl(aes(label=V3),
             method=list(cex=.6, 'last.points', 'filter', 'my.bumpup')) +
     scale_x_continuous(name='Negative log association p-value', breaks=cutoffs, labels=function(x) {sprintf("%.1f", x)}) +
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

args <- commandArgs(TRUE)
X <- read.table(args[1], header=FALSE, sep=' ')
w <- 4 * length(table(X$V2))
X$x <- -pchisq(X$V4, df=1, log.p=TRUE, lower.tail=FALSE)
cutoffs <- unique(X$V4)
if (length(args) > 1) {
  ## Plot enrichment of disjoint intervals between ticks
  offset <- 10 ^ (-.5 * log10(2))
} else {
  offset <- 1
}
margin <- 4
P <- 
panelsize <- 4
aspect=1
h <- panelsize * length(table(X$V1))
w <- (aspect * panelsize + margin) * length(table(X$V2))
pdf(file=gsub("in", "pdf", args[1]), height=4 * length(table(X$V1)), width=w)
t <- ggplot_gtable(ggplot_build(P))
t$layout$clip[t$layout$name == "panel"] <- "off"
grid.draw(t)
dev.off()
