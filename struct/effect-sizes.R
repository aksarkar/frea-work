library(directlabels)
library(ggplot2)
library(grid)
library(gridExtra)
library(scales)
library(zoo)
library(Cairo)

source("~/code/enr/plot/theme_nature.R")

# https://stackoverflow.com/a/16964861
facet_wrap_labeller <- function(gg.plot, labels=NULL) {
  g <- ggplotGrob(gg.plot)
  gg <- g$grobs
  strips <- grep("strip_t", names(gg))

  for(ii in seq_along(labels)) {
    modgrob <- getGrob(gg[[strips[ii]]], "strip.text",
                       grep=TRUE, global=TRUE)
    gg[[strips[ii]]]$children[[modgrob$name]] <- editGrob(modgrob,label=labels[ii])
  }

  g$grobs <- gg
  class(g) = c("arrange", "ggplot", class(g))
  g
}

density_by_bin <- function(x, aes='y=..density..') {
  (ggplot(data[!is.na(data$bin),], aes(x=odds_ratio, color=category, label=category)) +
   geom_density(aes_string(aes), size=I(.35 / ggplot2:::.pt)) +
   scale_y_continuous(name='Density', labels=comma) +
   scale_color_brewer(name='Category', palette='Dark2') +
   expand_limits(y=0) +
   facet_wrap(~ bin, ncol=4, scales='free_y') +
   theme_nature +
   theme(panel.margin=unit(4, 'mm'),
         plot.margin=unit(rep(1, 4), 'mm'),
         legend.position='bottom'))
}

effect_size <- function(x) {
  (density_by_bin(x, aes='y=..count..') +
   geom_vline(xintercept=1, linetype='dashed', size=I(.35 / ggplot2:::.pt)) +
   scale_x_continuous(name='Odds ratio', limits=c(0, 2)))
}

daf <- function(x) {
  (density_by_bin(x) +
   scale_x_continuous(name='Derived allele frequency', limits=c(0, 1)))
}

daf_by_effect_size_by_bin <- function(x) {
  (ggplot(data[!is.na(data$bin),], aes(x=daf, y=odds_ratio)) +
   geom_point(size=1, alpha=.25) +
   geom_hline(yintercept=1, linetype='dashed', size=I(.35 / ggplot2:::.pt)) +
   geom_vline(xintercept=.5, linetype='dashed', size=I(.35 / ggplot2:::.pt)) +
   scale_x_continuous(name='Derived allele frequency', limits=c(0, 1)) +
   scale_y_continuous(name='Odds ratio') +
   scale_fill_gradient(low='gray80', high='red') +
   facet_grid(bin ~ category) +
   theme_nature +
   theme(panel.margin=unit(4, 'mm'),
         plot.margin=unit(rep(1, 4), 'mm')))
}

args <- commandArgs(TRUE)
data <- read.table(args[1], col.names=c('category', 'daf', 'odds_ratio', 'logp'), row.names=NULL, sep=' ')
data$category <- factor(data$category, levels=c('non-synonymous', 'coding', 'ambiguous', 'non-coding'))
bins <- c(seq(0, 8, 2), 10 * seq(1, 3), max(data$logp))
data$bin <- cut(data$logp, breaks=bins, right=FALSE)
data$bin <- factor(data$bin, levels=rev(levels(data$bin)))
p <- daf_by_effect_size_by_bin(data)
Cairo(type='png', file=sub(".in$", ".png", args[1]), width=189, height=277, units='mm', dpi=300)
print(p)
dev.off()
