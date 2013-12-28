library(ggplot2)

args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=FALSE, sep=' ')
p <- (ggplot(d, aes(x=V3, y=V4, group=V3)) +
      geom_boxplot() +
      scale_x_continuous(name='Rank', breaks=seq(0, 300000, 50000)) +
      scale_y_continuous(name='Independent loci', breaks=seq(0, 150000, 30000)) +
      facet_grid(V1 ~ .) +
      theme_bw())
pdf()
print(p)
dev.off()

