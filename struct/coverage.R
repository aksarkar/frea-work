library(ggplot2)
library(scales)
library(Cairo)

source("~/code/enr/plot/color_roadmap.R")
source("~/code/enr/plot/theme_nature.R")

args <- commandArgs(TRUE)
data <- read.table(args[1], header=0, sep=' ')
data$x <- factor(data$V1, levels=enh_cluster_ordering)
p <- (ggplot(data, aes(x=x, y=V2, color=x)) +
      geom_point(size=I(1)) +
      labs(x='Cluster', y='Base pairs covered') +
      scale_y_continuous(labels=comma) +
      scale_color_roadmap +
      theme_nature +
      theme(axis.text.x=element_blank(),
            legend.position='none'))
Cairo(type='pdf', name=sub(".txt$", ".pdf", args[1]), dpi='auto', width=193, height=50, units='mm')
print(p)
dev.off()
