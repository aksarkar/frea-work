# Plot null distributions of rotations versus permutations
# Author: Abhishek Sarkar <aksarkar@mit.edu>

library(Cairo)
library(ggplot2)

args <- commandArgs(TRUE)
df <- melt(read.csv(args[1]))
p <- (qplot(data=df, x=value, xlab='Rank sum', ylab='Density', geom='density',
            color=factor(variable)) +
      scale_color_manual(name='Null distribution',
                         values=c('Rotations'='black', 'Permutations'='darkblue')) +
      scale_x_continuous(breaks=seq(floor(min(df$value) / 1e8) * 1e8,
                           ceiling(max(df$value) / 1e8) * 1e8, 1e8)) +
      geom_vline(color='red', xintercept=6718507428) +
      theme_bw())
CairoPDF(file='out.pdf', width=8, height=4.5)
print(p)
dev.off()
