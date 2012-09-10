library(ggplot2)
args <- commandArgs(TRUE)
d <- read.table(args[1], header=FALSE, sep=' ')
p <- (qplot(data=d, x=V5, y=V4, geom='hex',
           binwidth=c(1e4, .01),
           xlab='Average distance to closest TSS (bp)',
           ylab='Proportion of cell types active in') +
      scale_fill_gradient(low='#eeeeee', high='black') +
      theme_bw())
pdf()
print(p)
dev.off()
