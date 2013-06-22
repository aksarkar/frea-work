# RR plots for enrichment down the rank list
# Author: Abhishek Sarkar <aksarkar@mit.edu>
library(directlabels)
library(ggplot2)
library(grid)
library(plyr)
library(Cairo)

scale_state_aggregate <-
  scale_color_manual(name='State (aggregate)',
                     values=c('promoter' = '#ff0000', 'enhancer' = '#faca00',
                       'insulator' = '#09befe', 'transcribed' = '#00b050',
                       'repressed' = '#7f7f7f', 'other' = 'black',
                       'poly-A+-RNA-seq' = '#005c1f', 'diff-expressed' = '#0060a0'))

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

rrplot <- function(X, name, zero, cutoff=30000) {
  stopifnot(is.data.frame(X))
  stopifnot(is.numeric(zero))
  stopifnot(cutoff > 0)
  return(ggplot(X[X$total <= cutoff, ],
                aes(x=total, y=y, color=celltype)) +
         geom_line(size=I(.25)) +
         geom_hline(yintercept=zero, color="black", size=I(.25)) +
         geom_dl(aes(label=celltype),
                 method=list(cex=.8, "first.points", filter(10), "my.bumpup")) +
         facet_grid(phenotype ~ feature, scale="free") +
         scale_x_continuous(name="Rank", limits=c(-.6 * cutoff, cutoff),
                            breaks=seq(0, cutoff, cutoff / 4)) +
         scale_y_continuous(name=name) +
         theme_bw() +
         theme(strip.background=element_blank(),
               legend.position="none"))
}

args <- commandArgs(TRUE)
D <- read.csv(args[1], header=FALSE)
colnames(D) <- c("total", "phenotype", "celltype", "feature", "count", "expected")
if (length(args) > 1) {
  E <- read.csv(args[2], header=FALSE)
  colnames(E) <- c(colnames(D), "rep")
}
panelsize <- 15
h <- panelsize * length(table(D$phenotype))
w <- 1.6 * panelsize * length(table(D$feature))
Cairo(file=gsub(".in", ".pdf", args[1]), type="pdf", width=w, height=h, units="cm", dpi="auto")
rrplot(fold(D), "Fold enrichment", 1, 20000)
if (exists("E")) {
  rrplot(z.score(D, E, fold), "Z-score (fold enrichment)", 0, 150000)
}
warnings()
dev.off()
