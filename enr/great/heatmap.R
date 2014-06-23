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
    print(ggplot(X, aes(x=path, y=cell, fill=fold)) +
          geom_tile() +
          scale_fill_gradient(name='Fold enrichment', low='white', high='red') +
          scale_x_discrete(name=X$V4[1]) +
          scale_y_discrete(name='Cell type') +
          coord_fixed(ratio=1) +
          coord_flip() +
          theme_bw() +
          theme(legend.position='bottom',
                text=element_text(size=10),
                axis.text.x=element_text(angle=-90,hjust=0,vjust=.5),
                ## axis.text=element_blank(),
                axis.ticks=element_blank(),
                panel.grid=element_blank()))
}

args <- commandArgs(TRUE)
d <- read.delim(gzfile(args[1]), header=0)
load(args[2])
d$cell <- strtrim(d$V1, 40)
d$path <- strtrim(d$V2, 50)
d$fold <- d$V3
m <- cast(d, cell~path, value='fold', fun.aggregate=length)
d$cell <- factor(d$cell, levels=vapply(sample_info$name[order(sample_info$position)], function(x) {gsub(' ', '_', x)}, ""))
path.h <- hclust(dist(t(m)))
d$path <- factor(d$path, levels=path.h$labels[rev(path.h$order)])
d$clust <- cutree(path.h, k=64)[d$path]
limits=c(0, max(d$fold))

pdf(height=17, width=11)
## d_ply(d, .(clust), heatmap)
heatmap(d)
dev.off()
