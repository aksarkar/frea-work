library(ggplot2)
library(grid)
library(plyr)
library(reshape)
library(Cairo)

cluster.labels <- function(m) {
  h <- hclust(dist(m))
  h$labels[h$order]
}

heatmap <- function(X, limits=NULL) {
  if (is.null(limits)) {
    limits <- c(0, max(X$fold))
  }
  print(ggplot(X, aes(x=cell, y=path, fill=fold)) +
      geom_tile() +
      scale_fill_gradient(name='Log fold enrichment', low='white', high='red', limits=limits) +
      scale_x_discrete(name='Cluster', limits=cell.labels) +
      scale_y_discrete(name='GO category', limits=path.labels[path.labels %in% X$path]) +
      coord_equal() +
      theme_bw() +
      theme(text=element_text(size=6),
            axis.text.x=element_text(angle=90, hjust=1),
            axis.ticks=element_blank(),
            panel.grid=element_blank()))
}

args <- commandArgs(TRUE)
d <- read.delim(args[1], header=0)
d$cell <- strtrim(d$V1, 40)
d$path <- strtrim(d$V2, 50)
d$fold <- log10(d$V3)
m <- cast(d, cell~path, value='fold')
cell.labels <- cluster.labels(m)
## cell.labels <- cell.labels[cell.labels %in% c(17, 5, 43, 27, 19)]
path.h <- hclust(dist(t(m)))
path.labels <- path.h$labels[path.h$order]
d$clust <- cutree(path.h, k=16)[d$path]
limits=c(0, max(d$fold))
Cairo(type='svg', width=8.5, height=11, units='in')
## heatmap(d[d$clust %in% c(2,3,4,5,6,7,8,11,12,13,15,16),])  # enhancers
heatmap(d[d$clust %in% c(1,5,9,15),])
dev.off()
