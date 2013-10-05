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

hilbert.plot <- function(d, e) {
  n <- 16384
  to.uv <- function(x) {
    ldply(x$start, x.to.uv, n)
  }
  d <- ddply(d, .(chr, celltype), to.uv)
  e <- ddply(e, .(chr, score), to.uv)
  return(ggplot(d, aes(x=u, y=v)) +
         geom_point(aes(color=celltype), size=I(1), shape=I(3)) +
         geom_point(data=e, aes(alpha=score), size=I(1), shape=I(15)) +
         scale_x_continuous(limits=c(0, n), expand=c(0, 0)) +
         scale_y_continuous(limits=c(0, n), expand=c(0, 0)) +
         scale_color_brewer(name="Cell type", palette="Spectral",
                            guide=guide_legend(ncol=1, title.position="top")) +
         scale_alpha(name="SNP density", range=c(0, .8),
                     guide=guide_legend(ncol=1, title.position="top")) +
         coord_equal() +
         facet_wrap(~ chr, nrow=3) +
         theme_bw() + 
         theme(text=element_text(size=8),
               axis.text.x=element_text(angle=-90, hjust=0, vjust=0.5),
               axis.ticks=element_blank(),
               axis.title=element_blank(),
               legend.position='bottom',
               strip.background=element_blank()))
}

args <- commandArgs(TRUE)
d <- read.delim(args[1], header=FALSE)
colnames(d) <- c('chr', 'start', 'end', 'name', 'score', 'phenotype', 'feature', 'celltype')
e <- read.delim(args[2], header=FALSE)
colnames(e) <- c('chr', 'start', 'end', 'name', 'score')
d$chr <- factor(d$chr,
                levels=c(vapply(seq(1, 22), function(x) sprintf('chr%d', x), ''),
                  'chrX'))
p <- hilbert.plot(d, e)
svg(file=sub(".in", ".svg", args[1]), width=11, height=8.5)
print(p)
dev.off()
