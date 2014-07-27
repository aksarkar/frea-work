library(ggplot2)
library(Cairo)
args <- commandArgs(TRUE)
d <- read.table(args[1], header=0, sep=' ')
d$V3 <- d$V3 + rnorm(length(d$V3)) * 1e-7
p <- (ggplot(d, aes(x=V2, y=V3)) +
      geom_violin(aes(group=V2)) +
      stat_summary(fun.y=mean, geom='line', color='red') +
      geom_hline(yintercept=as.numeric(args[2]), size=I(.25)) +
      scale_y_continuous(name=expression(h[g]^2)) +
      scale_x_log10(name='Top n tags', labels=unique(d$V2), breaks=unique(d$V2)) +
      theme_bw() +
      theme(text=element_text(size=10),
            axis.text.x=element_text(angle=-90, hjust=0, vjust=.5),
            axis.title.y=element_text(angle=0),
            panel.grid.minor=element_blank(),
            strip.background=element_blank()))
Cairo(type='pdf', file=sub('.txt', '.pdf', args[1]), dpi='auto', width=160, height=100, units='mm')
print(p)
dev.off()
