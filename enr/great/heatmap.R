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
      scale_fill_gradient(name='Fold enrichment', low='white', high='red', limits=limits) +
      scale_x_discrete(name='Cell type', limits=cell.labels) +
      scale_y_discrete(name='GO biological process', limits=path.labels[path.labels %in% X$path]) +
      coord_equal() +
      theme_bw() +
      theme(text=element_text(size=6),
            axis.text.x=element_text(angle=-90, hjust=0, vjust=.5),
            axis.ticks=element_blank(),
            panel.grid=element_blank()))
}

args <- commandArgs(TRUE)
d <- read.delim(args[1], header=0)
d$cell <- strtrim(d$V1, 40)
d$path <- strtrim(d$V2, 50)
d$fold <- d$V3
m <- cast(d, cell~path, value='fold')
cell.labels <- cluster.labels(m)
path.h <- hclust(dist(t(m)))
path.labels <- path.h$labels[path.h$order]
d$clust <- cutree(path.h, k=16)[d$path]
limits=c(0, max(d$fold))

svg()
d_ply(d, .(clust), heatmap)
dev.off()
