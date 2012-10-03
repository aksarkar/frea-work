library(Cairo)
library(ggplot2)
args <- commandArgs(TRUE)
d <- read.csv(args[1])
p <- (qplot(data=d, x=cutoff, xlab='Top n SNPs', y=-log10(fisher),
            ylab="Negative log-transformed Fisher's exact test p-value",
            geom='line', size=I(.25)) +
      facet_grid(disease2 ~ disease1) +
      scale_x_continuous(breaks=seq(0, 20000, 2500)) +
      theme_bw() +
      opts(axis.text.x=theme_text(angle=-90, hjust=0),
           strip.background=theme_rect(fill=NA, colour=NA)))
CairoPDF(file='out.pdf', height=8, width=8)
print(p)
dev.off()
