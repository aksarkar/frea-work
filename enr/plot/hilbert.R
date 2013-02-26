library(bitops)
library(ggplot2)
library(grid)
library(plyr)

x.to.uv <- function(x, n) {
  u <- 0
  v <- 0
  s <- 1
  while (s < n) {
    rotate.x <- bitAnd(bitShiftR(x, 1), 1)
    rotate.y <- bitAnd(bitXor(x, rotate.x), 1)
    if (!rotate.y) {
      if (rotate.x) {
        u = s - 1 - u;
        v = s - 1 - v;
      }
      t <- v
      v <- u
      u <- t
    }
    u <- u + s * rotate.x
    v <- v + s * rotate.y
    x <- bitShiftR(x, 2)
    s <- bitShiftL(s, 1)
  }
  c(u=u, v=v)
}

hilbert.plot <- function(d) {
  n <- max(ddply(d, .(chr), function(x) {
    bitShiftL(1, ceiling(log(sqrt(max(x$start))) / log(2)))
  })$V1)
  to.uv <- function(x) {
    cbind(data.frame(score=x$score), ldply(x$start, x.to.uv, n))
  }
  d <- ddply(d, .(chr), to.uv)
  return(ggplot(d, aes(x=u, y=v, color=score)) +
         geom_point(size=I(1)) +
         scale_x_continuous(expand=c(0, 0)) +
         scale_y_continuous(expand=c(0, 0)) +
         scale_color_gradient(name='Beta', low='yellow', high='red') +
         expand_limits(x=0, y=0) +
         coord_equal() +
         facet_wrap(~ chr) +
         theme_bw() + 
         theme(text=element_text(size=8),
               axis.text.x=element_text(angle=90, hjust=1),
               axis.ticks=element_blank(),
               axis.title=element_blank(),
               strip.background=element_blank()))
}

args <- commandArgs(TRUE)
d <- read.delim(args[1], header=FALSE)
colnames(d) <- c('chr', 'start', 'end', 'name', 'score')
p <- hilbert.plot(d)
pdf(width=12, height=12)
print(p)
dev.off()
