library(ggplot2)
library(gtable)
library(grid)
library(plyr)
library(Cairo)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=FALSE, sep=' ')
p <- dlply(subset(d, V4 > 1 & V5 < 3e-5), .(V2),
           function(X) {
               Y <- X[order(X$V4) < 15,]
               return(ggplot(Y, aes(y=V3, x=log10(V4))) +
                      geom_point(size=I(1)) +
                      scale_x_continuous(name='Log of odds ratio', limits=c(0, 2),
                                         breaks=c(0, 1, 2), expand=c(0, 0)) +
                      scale_y_discrete(limits=Y[order(Y$V4),]$V3) +
                      facet_grid(V2 ~ .) +
                      coord_fixed(ratio=.35) +
                      theme_bw() +
                      theme(text=element_text(size=8),
                            rect=element_blank(),
                            axis.title=element_blank(),
                            panel.background=element_blank(),
                            panel.grid.minor=element_blank(),
                            plot.background=element_blank(),
                            strip.background=element_blank()))
           })
Cairo(type='pdf', dpi=96, height=300, width=178, units='mm')
grid.draw(do.call(rbind, c(llply(p, ggplotGrob), size='last')))
dev.off()
