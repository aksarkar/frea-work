## Visualize annotations in LD neighborhoods of GWAS markers
## Author: Abhishek Sarkar
## Arguments: INPUTFILE
## Input format: CSV (y, start, end, annotation) 
library(ggplot2)
library(Cairo)

scale_state_aggregate <-
  scale_fill_manual(name='State (aggregate)',
                     values=c('promoter' = '#ff0000', 'enhancer' = '#faca00',
                       'insulator' = '#09befe', 'transcribed' = '#00b050',
                       'repressed' = '#7f7f7f', 'other' = 'black',
                       'coding' = '#005c1f', 'protein_coding' = '#0060a0',
                       'known_lincRNA' = '#a000a0', 'putative_lincRNA' = '#ff00ff'))

args <- commandArgs(TRUE)
rects <- read.csv(args[1], header=0)
rects$V5 <- rects$V1 + 1
CairoPNG(width=2048, height=8192)
print(qplot(data=rects, geom='rect', xmin=V2, xmax=V3, ymin=V1, ymax=V5,
            fill=factor(V4)) +
      geom_vline(xintercept=0, color='black') +
      scale_state_aggregate +
      theme_bw() +
      opts(axis.text.y=theme_blank(), axis.ticks=theme_blank()))
dev.off()
