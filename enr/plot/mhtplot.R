# Manhattan plot of GWAS p-values
# Author: Abhishek Sarkar

library(ggplot2)
library(scales)
library(Cairo)

source("~/code/enr/plot/theme_nature.R")

args <- commandArgs(TRUE)
df <- read.table(args[1], header=FALSE, col.names=c('chr', 'pos', 'logp', 'category'))
lens <- c(0, 245203898, 243315028, 199411731, 191610523, 180967295, 170740541,
          158431299, 145908738, 134505819, 135480874, 134978784, 133464434,
          114151656, 105311216, 100114055, 89995999, 81691216, 77753510,
          63790860, 63644868, 46976537, 49476972, 152634166, 50961097)
xs <- cumsum(lens)
df$x = xs[df$chr] + df$pos
p <- (ggplot(df, aes(x=x, y=logp, color=factor(category))) +
      geom_point(size=I(.25), alpha=.5) +
      labs(x='Position', y=expression(-log(p))) +
      ## scale_x_continuous(breaks=xs, labels=seq(1, 25)) +
      scale_x_continuous(labels=comma) +
      scale_color_manual(values=c('1' = 'red', '0' = 'black')) +
      theme_nature)
Cairo(type='pdf', file=sub('.in$', '.pdf', args[1]), width=189, height=50, units='mm', dpi=96)
print(p)
dev.off()
