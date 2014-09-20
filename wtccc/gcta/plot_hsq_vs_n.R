library(ggplot2)
library(grid)
library(Cairo)

source("~/code/enr/plot/theme_nature.R")

args <- commandArgs(TRUE)
d <- read.table(args[1], header=0, sep=' ')
d$V3 <- d$V3 + rnorm(length(d$V3)) * 1e-7
p <- (ggplot(d, aes(x=V2, y=V3)) +
      geom_violin(aes(group=V2), size=I(.25)) +
      stat_summary(fun.y=mean, geom='line', color='red', size=I(.25)) +
      geom_hline(yintercept=as.numeric(args[2]), size=I(.25)) +
      scale_y_continuous(name=expression(h[g]^2)) +
      scale_x_log10(name='Top n tags', labels=unique(d$V2), breaks=unique(d$V2)) +
      theme_nature +
      theme(axis.title.y=element_text(angle=0),
            axis.text=element_text(size=5),
            axis.text.x=element_text(angle=-90, hjust=0, vjust=.5),
            plot.margin=unit(c(1, 0, 0, 0), 'mm')
            ))
Cairo(dpi='auto', file=sub('.txt', '.pdf', args[1]), height=37, type='pdf', units='mm', width=60)
print(p)
dev.off()
