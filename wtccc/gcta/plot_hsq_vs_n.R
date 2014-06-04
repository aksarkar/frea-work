library(ggplot2)
library(Cairo)
args <- commandArgs(TRUE)
d <- read.table(args[1], header=0, sep=' ')
h2g <- data.frame(V5=c('t1d', 'permuted', 'random'), V3=c(.336, 0, .336))
d$V3 <- d$V3 + rnorm(length(d$V3)) * 1e-7
p <- (ggplot(d, aes(x=V2, y=V3)) +
      geom_violin(aes(group=V2)) +
      stat_summary(fun.y=mean, geom='line', color='red') +
      geom_hline(data=h2g, size=I(.25)) +
      scale_y_continuous(name=expression(h[g]^2)) +
      scale_x_log10(name="Top n tags", labels=unique(d$V2), breaks=unique(d$V2)) +
      theme_bw() +
      theme(text=element_text(size=10),
            axis.text.x=element_text(angle=-90, hjust=0, vjust=.5),
            axis.title.y=element_text(angle=0),
            panel.grid.minor=element_blank(),
            strip.background=element_blank()))
Cairo(type="svg", file=sub(".txt", ".svg", args[1]), dpi="auto", width=160, height=100, units='mm')
print(p)
dev.off()
Cairo(type="png", file=sub(".txt", ".png", args[1]), dpi=200, width=160, height=100, units='mm')
print(p)
dev.off()
