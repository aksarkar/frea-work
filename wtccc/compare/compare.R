# Compare WTCCC1 imputed p-values with our LD-expanded p-values
# Author: Abhishek Sarkar <aksarkar@mit.edu>
library(Cairo)
library(grid)
library(ggplot2)
library(hexbin)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=TRUE, sep=' ')
d$diff = d$expanded - d$imputed
p1 <- (qplot(data=d, x=diff, color=factor(study),
             xlab='Representative p-value minus imputed p-value',
             geom='freqpoly', binwidth=.05) +
       scale_color_brewer(palette='Set3', name='Study') +
       theme_bw())
p2 <- hexbin(d$imputed, d$expanded, xlab='Imputed p-value',
             ylab='Representative p-value', xbins=100)
pdf(file='compare-dist.pdf', height=6, width=6)
print(p1)
dev.off()
pdf(file='compare-hexbin.pdf', height=6, width=6)
plot(p2)
dev.off()
