## Visualize distribution of functional elements around SNPs
## Author: Abhishek Sarkar

library(Cairo)
library(ggplot2)

args <- commandArgs(TRUE)
d <- read.table(gzfile(args[1]), header=0, sep=' ')
d$y = rep(0, len=length(d$V1))
p <- (ggplot(d, aes(x=as.numeric(factor(V3)), y=y, ymin=V1, ymax=V2)) +
      geom_linerange(position=position_jitter(width=.1)) +
      scale_x_continuous(name='Cell type', breaks=as.numeric(factor(d$V3)),
                         labels=factor(d$V3)) +
      scale_y_continuous(name='Relative position') +
      coord_flip() +
      theme_bw())
Cairo(file=sub('.in.gz$', '.png', args[1]), type='png', dpi=96)
print(p)
dev.off()
