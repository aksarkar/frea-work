library(ggplot2)
library(Cairo)
args <- commandArgs(TRUE)
d <- read.table(args[1], header=0, sep=' ')
p <- (ggplot(d, aes(x=V2, y=V3, group=V2)) +
      geom_boxplot(width=.1, outlier.size=1) +
      scale_y_continuous(name="Pruned samples") +
      scale_x_log10(name="Top n tags", labels=unique(d$V2), breaks=unique(d$V2)) +
      theme_bw() +
      theme(text=element_text(size=10),
            axis.text.x=element_text(angle=-90, hjust=0, vjust=.5),
            panel.grid.minor=element_blank(),
            strip.background=element_blank()))
Cairo(type="pdf", file=sub(".txt", ".pdf", args[1]), dpi=96, width=160, height=100, units='mm')
print(p)
dev.off()
