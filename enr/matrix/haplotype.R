library(ggplot2)
library(scales)
library(Cairo)

source("~/code/enr/plot/theme_nature.R")

args <- commandArgs(TRUE)
data <- read.table(gzfile(args[1]), header=FALSE, sep=' ')
data$x <- with(data, pmin(V2, V3))
data$y <- with(data, pmax(V2, V3))
p <- (ggplot(data, aes(x=x, y=y, color=V5)) +
      labs(x='Position', y='Position', title=sub('.txt.gz$', '', args[1])) +
      geom_point(size=I(.25)) +
      scale_x_continuous(labels=comma, expand=c(0, 0)) +
      scale_y_continuous(labels=comma, expand=c(0, 0)) +
      ## scale_color_gradient(limits=c(.1, 1), low='#ffeda0', high='#f03b20') +
      scale_color_gradient2(low=muted('blue'), mid='white', high=muted('red')) +
      coord_fixed() +
      theme_nature +
      theme(plot.title=element_text(size=7)))
Cairo(type='png', height=89, width=89, units='mm', dpi=300, bg='white')
print(p)
dev.off()
