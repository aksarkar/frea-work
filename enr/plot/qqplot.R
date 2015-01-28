library(directlabels)
library(ggplot2)
library(plyr)
library(Cairo)

source('~/code/enr/plot/color_roadmap.R')
source('~/code/enr/plot/theme_nature.R')

qqplot <- function(data, theoretical) {
  # Expected QQ line under the null
  quantiles <- c(.25, .75)
  y <- quantile(data$observed, quantiles)
  x <- theoretical(quantiles)
  m <- diff(y) / diff(x)
  b <- data$observed[1L] - m * data$expected[1L]
  m <- 1
  b <- 0
  (ggplot(data, aes(x=expected, y=observed, color=category)) +
   geom_line(aes(x=all_expected, group='all'), size=I(.5 / ggplot2:::.pt), color='black') +
   geom_line(size=I(.5 / ggplot2:::.pt)) +
   geom_abline(slope=m, yintercept=b, linetype='dashed', size=I(.35 / ggplot2:::.pt)) +
   geom_dl(aes(label=category), method=list(cex=5/16, 'last.qp')) +
   scale_color_brewer(palette='Dark2') +
   labs(x="Expected -log10(p)", y="Observed -log10(p)") +
   coord_cartesian(xlim=c(0, 5), ylim=c(0, 10)) +
   theme_nature +
   theme(plot.margin=unit(rep(3, 4), 'mm')))
}

args <- commandArgs(TRUE)
data <- read.table(args[1], col.names=c('category', 'observed'), row.names=NULL)
data <- data[order(data$observed, decreasing=TRUE),]
data$all_expected <- -log10(ppoints(data$observed))
subsampled <- ddply(data, .(category), function(x) {x[seq(1, nrow(x), nrow(x) / 3000),]})
with_expected <- ddply(subsampled, .(category), transform, expected=-log10(ppoints(observed)))
Cairo(type='pdf', file=sub(".in.gz$", ".pdf", args[1]), dpi='auto', width=83, height=83, units='mm')
qqplot(with_expected, function(x) {x})
dev.off()
