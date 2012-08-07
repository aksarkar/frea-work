library(Cairo)
library(ggplot2)
library(ggdendro)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=FALSE)
h <- hclust(dist(as.matrix(cast(d, V1 ~ V2, value='V3'))))
labels <- unique(d$V1)[h$order]

q <- ggdendrogram(h, leaf_labels=FALSE, rotate=TRUE)
CairoPDF(file='tree.pdf', height=16, width=2)
print(q)
dev.off()

p <- (qplot(data=d, x=V1, y=V2, fill=V3, geom='tile',
            xlab='Cell type (hierachically clustered)',
            ylab='Cell type (hierarchically clustered)') +
      scale_x_discrete(limits=labels) +
      scale_y_discrete(limits=labels) +
      scale_fill_gradient(low='white', high='red', name='bp in overlapping DHS regions') +
      coord_equal() +
      theme_bw() +
      opts(axis.text.x=theme_text(angle=-90, hjust=0)))
CairoPDF(file='heatmap.pdf', height=16, width=16)
print(p)
dev.off()
