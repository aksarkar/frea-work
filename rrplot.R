library(ggplot2)
library(Cairo)
args <- commandArgs(TRUE)
D <- read.csv(args[1])
D$denom = rep(tail(D$expected, length(table(D$feature))), length.out=length(D$count))
D$dev <- (D$count - D$expected) / D$denom

scale_cell_type <-
  scale_color_manual(name='Cell type',
                     values=c('H1' = '#3e8f54', 'Huvec' = '#e2560f',
                       'GM12878' = '#9e320e', 'HepG2' = '#1182bb',
                       'HMEC' = '#3fe7d2', 'K562' = '#3c11bc',
                       'NHEK' = '#bb1d94', 'HSMM' = '#e0c31b',
                       'NHLF' = '#806cd3'))

scale_state_aggregate <-
  scale_color_manual(name='State (aggregate)',
                     values=c('promoter' = '#ff0000', 'enhancer' = '#faca00',
                       'insulator' = '#09befe', 'transcribed' = '#00b050',
                       'repressed' = '#7f7f7f', 'other' = 'black'))

p <- (qplot(data=subset(D, total < 150000), x=total, y=dev, geom='line',
            size=I(.5), color=feature, xlab='Rank', ylab='Deviation') + 
      geom_hline(yintercept=0) + 
      scale_state_aggregate +
      scale_x_continuous(breaks=c(0, 50000, 100000, 150000)) +
      theme_bw())
Cairo(type='pdf', file='out.pdf', dpi=96)
print(p)
dev.off()
