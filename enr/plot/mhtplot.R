# Manhattan plot of GWAS p-values
# Author: Abhishek Sarkar

library(Cairo)
library(ggplot2)

scale_state_aggregate <-
  scale_color_manual(name='State (aggregate)',
                     values=c('promoter' = '#ff0000', 'enhancer' = '#faca00',
                       'insulator' = '#09befe', 'transcribed' = '#00b050',
                       'repressed' = '#7f7f7f', 'other' = 'black',
                       'coding' = '#005c1f', 'diffexpr' = '#0060a0'))

CairoFonts(regular='Dejavu Sans')
args <- commandArgs(TRUE)
df <- read.csv(args[1])
lens <- c(0, 245203898, 243315028, 199411731, 191610523, 180967295, 170740541,
          158431299, 145908738, 134505819, 135480874, 134978784, 133464434,
          114151656, 105311216, 100114055, 89995999, 81691216, 77753510,
          63790860, 63644868, 46976537, 49476972, 152634166, 50961097)
xs <- cumsum(lens)
df$x = xs[df$chr] + df$pos
p <- (qplot(data=df, x=x, xlab='Genomic position', y=-log10(p), geom='point',
            color=state, size=I(1)) +
      scale_x_continuous(breaks=xs, labels=seq(1, 25)) +
      scale_state_aggregate +
      geom_hline(color='red', yintercept=6) +
      scale_y_continuous(limits=c(0, ceiling(max(-log10(df$p))))) +
      opts(legend.position='none') +
      theme_bw())
postscript()
print(p)
dev.off()
