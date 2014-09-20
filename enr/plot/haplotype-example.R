library(ggplot2)
library(plyr)
library(Cairo)

source("~/code/enr/plot/theme_nature.R")

args <- commandArgs(TRUE)
data <- read.delim(gzfile(args[1]), header=FALSE)
relative.positions <- ddply(data, .variables=.(V4), .drop=FALSE,
                            .fun=function(X) {transform(X, pos=X$V2 - mean(X$V2))})
p <- (ggplot(relative.positions, aes(x=pos, y=factor(V4), color=factor(V6))) +
      geom_point(size=I(1)) +
      scale_color_manual(values=c('0' = 'gray70', '1' = 'red')) +
      scale_y_discrete(name='Tag SNP') +
      theme_nature +
      theme(axis.line=element_blank(),
            axis.ticks=element_blank(),
            axis.text.x=element_blank(),
            axis.title.x=element_blank(),
            legend.position='none',
            panel.grid.major=element_line(size=.25 / ggplot2:::.pt),
            panel.grid.major.x=element_blank()))
Cairo(type='pdf', dpi='auto', width=183, height=50, units='mm')
print(p)
dev.off()
