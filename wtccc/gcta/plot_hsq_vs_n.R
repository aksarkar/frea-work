library(ggplot2)
library(Cairo)
args <- commandArgs(TRUE)
d <- read.table(args[1], header=0, sep=' ')
d$V5 <- factor(d$V5, levels=c('T1D', 'Permuted'))
p <- (ggplot(d, aes(x=V2, y=V3)) +
      geom_pointrange(aes(ymin=V3 - V4, ymax=V3 + V4), size=I(.25),
                      position=position_jitter(width=.05)) +
      stat_smooth(method='lm', alpha=.25, color='#cc0000', se=0, size=.25) +
      geom_hline(yintercept=.33, size=I(.25)) +
      scale_y_continuous(name=expression(h[g]^2)) +
      scale_x_log10(name="Top n tags", labels=unique(d$V2), breaks=unique(d$V2)) +
      facet_grid(V5 ~ ., scales='free_y') +
      theme_bw() +
      theme(text=element_text(size=10),
            axis.text.x=element_text(angle=-90, hjust=0, vjust=.5),
            axis.title.y=element_text(angle=0),
            panel.grid.minor=element_blank(),
            strip.background=element_blank()))
Cairo(type="svg", file=sub(".txt", ".svg", args[1]), dpi=96, width=160, height=100, units='mm')
print(p)
dev.off()
