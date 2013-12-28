library(Cairo)
library(ggplot2)
library(grid)
library(gridExtra)
library(plyr)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=FALSE, sep=' ')
d$z <- (d$V4 - d$V5) / d$V6
p <- dlply(subset(d, z > 0), .(V1),
           function(X) {
               Y <- X[order(X$z, decreasing=TRUE) < 15,]
               ggplotGrob(ggplot(Y, aes(y=V2, x=z)) +
                      geom_point(size=I(1)) +
                      scale_x_continuous(name='z-score', limits=c(0, ceiling(max(Y$z)))) +
                      scale_y_discrete(limits=rev(Y[order(Y$z, decreasing=TRUE),]$V2)) +
                      facet_grid(V1 ~ .) +
                      coord_fixed(ratio=.1) +
                      theme_bw() +
                      theme(text=element_text(size=6),
                            rect=element_blank(),
                            axis.title=element_blank(),
                            panel.background=element_blank(),
                            panel.grid.minor=element_blank(),
                            plot.background=element_blank(),
                            strip.background=element_blank()))
           })
grid.draw(do.call(cbind, c(p, size='last')))
ggsave(file=sub('.in', '.pdf', args[1]))
