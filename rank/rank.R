library(Cairo)
library(ggplot2)
args <- commandArgs(TRUE)
d <- read.csv(args[1])
d$fold <- d$tt / d$exp
head(d)
p <- (qplot(data=d, x=cutoff, xlab='Top n SNPs',
            y=log10(fold), ylab='Log-transformed fold enrichment', geom='line', size=I(.25)) +
      facet_grid(disease2 ~ disease1) +
      ## scale_x_continuous(limits=c(0, 10000), breaks=seq(0, 10000, 2500)) +
      ## scale_y_continuous(limits=c(0, 10000), breaks=seq(0, 10000, 2500)) +
      theme_bw() +
      opts(axis.text.x=theme_text(angle=-90, hjust=0),
           strip.background=theme_rect(fill=NA, colour=NA)))
CairoPDF(file='out.pdf', height=8, width=8)
print(p)
dev.off()
