library(ggplot2)
library(Cairo)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=0)
p <- (ggplot(d, aes(x=V1, y=V2, ymin=(V2 - V3), ymax=(V2 + V3))) +
      scale_y_continuous(name="Proportion of heritability explained") +
      scale_x_log10(name="Top n tags", labels=d$V1, breaks=d$V1) +
      geom_line() +
      geom_errorbar(width=.05) +
      theme_bw() +
      theme(axis.text.x=element_text(angle=-90,hjust=0,vjust=.5),
            panel.grid.minor=element_blank()))
Cairo(type="pdf", name=sub(".in", ".pdf", args[1]), dpi=96)
print(p)
dev.off()
