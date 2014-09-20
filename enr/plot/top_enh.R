library(ggplot2)
library(grid)
library(reshape)
library(scales)
library(Cairo)

source('~/code/enr/plot/theme_nature.R')

args <- commandArgs(TRUE)
d <- read.delim(args[1], header=FALSE)
p <- (ggplot(d, aes(x=V1, y=V2, label=V3, fill=V3)) +
      geom_tile() +
      geom_text(size=(5 / ggplot2:::.pt)) +
      coord_equal() +
      scale_fill_gradient(low='white', high='red', guide='none') +
      scale_y_discrete(limits=rev(levels(d$V1))) +
      theme_nature +
      theme(text=element_text(size=5),
            axis.title=element_blank(),
            axis.text.x=element_text(angle=45, hjust=1, vjust=1)))
Cairo(type='pdf', file=sub('.txt$', '.pdf', args[1]), dpi='auto', height=89, width=89, units='mm')
print(p)
dev.off()
