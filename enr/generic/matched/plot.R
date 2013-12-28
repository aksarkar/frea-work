library(ggplot2)
library(grid)
library(plyr)
library(Cairo)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=FALSE, sep=' ')
P <- (ggplot(d, aes(x=V2, y=V3, fill=V4)) +
      geom_tile() +
      scale_x_discrete(name='Cell type', expand=c(0, 0)) +
      scale_y_discrete(name='Variant set', expand=c(0, 0)) +
      scale_fill_gradient2(name='Odds ratio', low='#0000cc', mid='white',
                           high='#cc0000', midpoint=1) +
      facet_wrap(~ V1, ncol=1) +
      coord_fixed(ratio=1) +
      theme_bw() +
      theme(strip.background=element_blank(),
            axis.text.x=element_text(angle=-90, hjust=0, vjust=.5)))
e <- d[d$V5 < (.05 / length(d$V5)),]
if (nrow(e) > 0) {
    P <- P + geom_text(data=e, label='*')
}
ggsave(file=sub('.in$', '.eps', args[1]))
