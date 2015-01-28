library(directlabels)
library(ggplot2)
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

dev <- function(X) {
  ddply(X, .(phenotype, feature, celltype), transform,
        y=(count - expected) / max(count))
}

normalize <- function(X) {
  ddply(X, .(phenotype), transform,
        z=y / max(y))

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

slice.plot <- function(X) {
  Y <- normalize(X[X$total == 10000, ])
  print(ggplot(Y, aes(x=phenotype, y=z, color=factor(celltype))) +
        geom_point(size=I(.5)) +          
        geom_hline(yintercept=0, color='black', size=I(.5 / ggplot2:::.pt)) +
        geom_dl(aes(label=celltype),
                method=list(cex=7/16, dl.trans(x=x+.5), 'bumpdown', top(10))) +
        scale_x_discrete(name='Feature') +
        scale_y_continuous(name='Relative cumulative deviation') +
        scale_color_roadmap +
        facet_wrap(~ phenotype, scales='free_x') +
        theme_nature +
        theme(legend.position='none',
              axis.text.x=element_blank(),
              axis.title.x=element_blank()))
}

args <- commandArgs(TRUE)
X <- read.csv(gzfile(args[1]), header=FALSE)
colnames(X) <- c('total', 'phenotype', 'celltype', 'feature', 'count', 'expected')
Cairo(type='pdf', width=89, height=60, units='mm', dpi='auto', family='Helvetica')
slice.plot(dev(X))
dev.off()
