library(ggplot2)
library(gtable)
library(grid)
library(plyr)
library(Cairo)

source('~/code/enr/plot/theme_nature.R')

my.panel <- function(X) {
  Y <- head(X[order(X$z, decreasing=TRUE),], n=15)
  ggplotGrob(ggplot(Y, aes(y=V3, x=z)) +
             geom_point(size=I(.5)) +
             scale_y_discrete(name='Cell type', limits=rev(Y$V3)) +
             scale_x_continuous(name='Enrichment z-score', limits=c(0, 12), breaks=seq(0, 12, 3)) +
             coord_fixed(ratio=1.618) +
             theme_nature +
             theme(text=element_text(size=5),
                   axis.title.y=element_blank()))
}

args <- commandArgs(TRUE)
data <- read.delim(args[1], header=FALSE)
data$z <- (data$V5 - data$V6) / sqrt(data$V7)
data$V2 <- factor(data$V2, levels=c('dhs', 'dgf', 'Enh'))
p <- dlply(data, .(V2), my.panel)
Cairo(type='pdf', dpi=96, height=100, width=63, units='mm')
grid.draw(do.call(cbind, c(p, size='last')))
dev.off()
