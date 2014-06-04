library(ggplot2)

args <- commandArgs(TRUE)
assoc <- read.delim(args[1])
y <- sort(assoc$WaldStat)
x <- qchisq(ppoints(y), df=1)
p <- (qplot(x, y, xlab='Expected', ylab='Observed', geom='hex') +
      scale_fill_gradient(low='white', high='red') +
      stat_abline(slope=1, yintercept=0) +
      coord_equal())
pdf()
print(p)
dev.off()
