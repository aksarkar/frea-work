library(Cairo)
library(ggplot2)

args <- commandArgs(TRUE)
D <- read.table(args[1], header=FALSE, sep=' ')
p <- (qplot(data=D, x=V2, y=V3, color=factor(V1), xlab='P-value (binned)',
            ylab='Average number of SNPs in LD (R^2 > .2)', geom='line',
            size=I(.5)) +
      scale_color_brewer(palette='Set2', name='Phenotype') +
      theme_bw())
CairoPDF(file='out.pdf')
print(p)
dev.off()
