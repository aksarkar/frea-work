library(ggplot2)
library(grid)
library(plyr)
library(Cairo)

source("~/code/enr/plot/theme_nature.R")

args <- commandArgs(TRUE)
d <- read.table(args[1], header=0, sep=' ')
total_expectation_variance <-
  ddply(d, .(V2), function (x) {
    expected_value <- mean(x$V3)
    total_variance <- var(x$V3) + mean(x$V4)
    data.frame(y=expected_value, ymin=expected_value - sqrt(total_variance),
               ymax=expected_value + sqrt(total_variance), m=median(x$V3))
  })
p <- (ggplot(total_expectation_variance, aes(x=V2, y=y, ymin=ymin, ymax=ymax)) +
      geom_pointrange(size=I(.35 / ggplot2:::.pt)) +
      geom_hline(yintercept=0, size=I(.25)) +
      geom_hline(yintercept=as.numeric(args[2]), size=I(.25), linetype='dashed') +
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
