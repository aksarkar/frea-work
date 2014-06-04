library(ggplot2)
library(grid)
library(reshape)
library(Cairo)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=FALSE, sep=' ')
h <- hclust(dist(as.matrix(cast(d, V1 ~ V2, value='V3'))))
labels <- unique(d$V1)[h$order]
p <- (ggplot(d, aes(x=V1, y=V2, label=V3, fill=V3)) +
      geom_tile() +
      geom_text(size=2.5) +
      coord_equal() +
      scale_fill_gradient(low='white', high='red', guide='none') +
      scale_x_discrete(name='Reference epigenome', limits=labels) +
      scale_y_discrete(name='Reference epigenome', limits=rev(labels)) +
      theme_bw() +
      theme(text=element_text(size=8),
            axis.text.x=element_text(angle=-45, hjust=0, vjust=1),
            plot.margin=unit(c(0,1,0,0), "in")))
Cairo(type='pdf', file=sub('.txt$', '.pdf', args[1]), dpi=96, height=6, width=8, units='in')
print(p)
dev.off()
