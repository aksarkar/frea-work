library(ggplot2)
library(gtable)
library(grid)
library(plyr)
library(Cairo)

args <- commandArgs(TRUE)
d <- read.delim(args[1], header=FALSE)
p <- dlply(d, .(V1),
           function(X) {
               ggplotGrob(ggplot(X, aes(y=V2, x=V3)) +
                          geom_point(size=I(1)) +
                          scale_x_log10(name='Log fold enrichment', limits=c(1, 50000), breaks=c(1, 10, 100, 1000, 10000), expand=c(0, 0)) +
                          scale_y_discrete(limits=X[order(X$V3),]$V2) +
                          facet_grid(V1 ~ .) +
                          coord_fixed(ratio=.35) +
                          theme_bw() +
                          theme(text=element_text(size=6),
                                rect=element_blank(),
                                axis.title=element_blank(),
                                axis.text.x=element_text(angle=-90, hjust=0, vjust=.5),
                                panel.background=element_blank(),
                                panel.grid.minor=element_blank(),
                                plot.background=element_blank(),
                                strip.background=element_blank(),
                                strip.text.y=element_text(angle=0, hjust=0, vjust=.5)))
           })
Cairo(file=sub(".hits.txt", ".svg", args[1]), type='svg', dpi=96, height=17, width=11, units='in')
grid.draw(do.call(rbind, c(p, size='last')))
dev.off()
