## Visualize distribution of functional elements around tagged markers
## Author: Abhishek Sarkar

library(Cairo)
library(ggplot2)

scale_state_aggregate <-
  scale_color_manual(name='State (aggregate)',
                     values=c('promoter' = '#ff0000', 'enhancer' = '#faca00',
                       'insulator' = '#09befe', 'transcribed' = '#00b050',
                       'repressed' = '#7f7f7f', 'other' = 'black',
                       'coding' = '#005c1f', 'protein_coding' = '#0060a0',
                       'known_lincRNA' = '#a000a0', 'putative_lincRNA' = '#ff00ff'))

args <- commandArgs(TRUE)
d <- read.csv(args[1], header=0)
CairoPDF(file='out.pdf', width=11, height=8.5)
print(qplot(data=d, x=V1, y=V3, geom='freqpoly', stat='identity',
            color=factor(V2), xlab='Relative position',
            ylab='Log-transformed count') +
      scale_state_aggregate +
      theme_bw())
dev.off()
