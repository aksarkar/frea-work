# Plot empirical PDF of SNPs in enhancers versus not
# Author: Abhishek Sarkar

library(Cairo)
library(ggplot2)

CairoFonts(regular='Dejavu Sans', symbol='Dejavu Sans')
args <- commandArgs(TRUE)
df <- read.delim(args[1], header=0)
df$enh = factor(df$V6, levels=c('0', '1'), labels=c('False', 'True'))
p <- (qplot(data=df, x=V5, xlab='p-value', ylab='density', geom='density',
           color=factor(enh)) +
      scale_color_manual(name='Enhancer?',
                         values=c('False'='grey', 'True'='black')) +
      theme_bw())
CairoPDF(file='out.pdf')
print(p)
dev.off()
