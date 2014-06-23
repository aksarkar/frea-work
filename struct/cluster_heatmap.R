library(ggplot2)
library(reshape2)

load('dat_summary.RData')
d <- subset(melt(dat_summary), value > .5)
p <- (ggplot(d, aes(x=factor(Var1), y=factor(Var2), fill=value)) +
      scale_x_discrete() +
      scale_y_discrete() +
      geom_tile() +
      scale_fill_gradient(low='white', high='red', guide='none') +
      coord_equal() +
      theme_bw() +
      theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=.5)))
pdf(height=17, width=11)
print(p)
dev.off()
