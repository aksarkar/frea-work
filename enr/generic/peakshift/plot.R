library(Cairo)
library(ggplot2)
library(reshape2)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=FALSE, sep=' ')
p <- (ggplot(d) +
      geom_point(aes(x=V2, y=V3/V4, color=(V3 / 1e4 < .01/50))) +
      scale_x_discrete(name='Functional category') +
      scale_y_continuous(name='Fold enrichment (peak shifting null)') +
      scale_color_manual(values=c('red', 'black'), guide='none') +
      theme_bw() +
      facet_grid(V1 ~ .) +
      theme(axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5),
            strip.background=element_blank()))
Cairo(name=sub('.in', '.pdf', args[1]), type='pdf', width=200, height=240, dpi=96, units='mm')
print(p)
dev.off()
