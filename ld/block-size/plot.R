library(Cairo)
library(ggplot2)

args <- commandArgs(TRUE)
D <- read.table(args[1], header=FALSE, sep=' ')
p <- (qplot(data=D, x=V2, y=V3, color=factor(V1), xlab='p-value',
            ylab='average number of SNPs in LD (R^2 > .2)', geom='line') +
      scale_color_brewer(palette='Set2') +
      theme_bw())
CairoPDF(file='out.pdf')
print(p)
dev.off()
