# Plot empirical PDF of SNPs in enhancers versus not
# Author: Abhishek Sarkar

library(Cairo)
library(ggplot2)

args <- commandArgs(TRUE)
df <- read.delim(args[1], header=0)
df$enh = factor(df$V6, levels=c('0', '1'), labels=c('False', 'True'))
p <- (qplot(data=df, x=V5, xlab='p-value', ylab='density', geom='density',
           color=factor(enh)) +
      scale_color_manual(values=c('False'='black', 'True'='#faca00'),
                         legend=FALSE) +
      theme_bw())
CairoPDF(file='out.pdf')
print(p)
dev.off()
