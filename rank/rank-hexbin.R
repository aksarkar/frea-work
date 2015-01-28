library(ggplot2)
library(grid)
library(Cairo)

source("~/code/enr/plot/theme_nature.R")

rank_hexbin <- function(data) {
  (ggplot(data[data$V4 < Inf & data$V5 < Inf,], aes(x=V4, y=V5)) +
   geom_hex() +
   geom_abline(slope=1, yintercept=0, linetype='dashed', size=.35 / ggplot2:::.pt) +
   scale_fill_gradient(name='Count', low='#fee8c8', high='#e34a33', trans='log', breaks=10 ** c(1, 3, 5), guide=guide_colorbar(barwidth=unit(30, 'mm'))) +
   ## xlim(0, 20)+
   ## ylim(0, 20) +
   theme_nature +
   theme(legend.position='bottom'))
}

args <- commandArgs(TRUE)
data <- read.delim(gzfile(args[1]), header=0)
labeller <- function(x) {paste(x, expression(-log(p)))}
p <- rank_hexbin(data) + labs(x=labeller(args[2]), y=labeller(args[3]))
pdf(file=sub('.bed.gz$', '.pdf', args[1]), width=83 / 25.4, height=88 / 25.4)
print(p)
dev.off()
