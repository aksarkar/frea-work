# Plot observed versus expected waiting times between genomic elements
# Author: Abhishek Sarkar

library(Cairo)
library(ggplot2)

CairoFonts(regular='Dejavu Sans')
args <- commandArgs(TRUE)
df <- read.delim(args[1], header=0)
df$V2 = rexp(length(df$V1), length(df$V1) / 3e9)
colnames(df) <- c('Observed', 'Expected')
d <- subset(melt(df), value >= 1)
p <- (qplot(data=d, x=log10(value), xlab='Log-transformed waiting time',
            y=log10(..count..), ylab='Log-transformed count',
            color=factor(variable),
            stat='bin', binwidth=.1, geom='line') +
      scale_color_manual(values=c('Observed'='red', 'Expected'='black'), legend=FALSE) +
      theme_bw())
CairoPDF(file='out.pdf', width=8, height=4.5)
print(p)
dev.off()
