library(ggplot2)
library(plyr)

args <- commandArgs(TRUE)
d <- read.delim(args[1], sep=' ', header=FALSE)
bg <- read.delim(args[2], sep=' ', header=FALSE)
e <- ddply(d, .(V1, V2, V3), function(x) x[which.max(x$V6),])
p <- (ggplot(e, aes(x=V2, y=V3, label=V4, fill=log(V6))) +
      scale_x_discrete(name='Feature', limits=c('dhs', 'distal_dhs', 'enh',
      'distal_enh', 'dhs+enh', 'dgf', 'distal_dgf', 'faire-seq',
      'distal_faire-seq')) +
      scale_y_discrete(name='Cell type') +
      scale_fill_gradient(name='Max log fold enrichment', low='white', high='red') +
      geom_tile(color='black') +
      coord_equal() +
      facet_grid(. ~ V1) +
      theme_bw() +
      theme(axis.text.x=element_text(size=6, angle=90, hjust=1),
            axis.text.y=element_text(size=6),
            panel.grid.major=element_blank(),
            strip.text.x=element_text(size=6),
            strip.background=element_blank())
      )
pdf(width=10, height=20)
print(p)
dev.off()
