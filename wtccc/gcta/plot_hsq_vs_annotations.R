library(ggplot2)
library(gridExtra)
library(Cairo)

source("~/code/enr/plot/color_roadmap.R")
source("~/code/enr/plot/theme_nature.R")

## https://stackoverflow.com/a/24241954
fancy_scientific <- function(l) {
    ## turn in to character string in scientific notation
    l <- format(l, scientific = TRUE)
    l <- gsub("0e\\+00","0",l)
    ## quote the part before the exponent to keep all the digits
    l <- gsub("^(.*)e", "'\\1'e", l)
    ## turn the 'e+' into plotmath format
    l <- gsub("e", "%*%10^", l)
    ## return this as an expression
    parse(text=l)
}

args <- commandArgs(TRUE)
d <- read.table(args[1], header=0)
d$V2 <- factor(d$V2, levels=enh_cluster_ordering)
p <- (ggplot(d, aes(x=V2, y=(V3), fill=V2)) +
      geom_bar(stat='identity') +
      scale_x_discrete(name="Enhancer module", drop=FALSE) +
      scale_y_continuous(name="Proportion of variance explained", breaks=seq(0, .04, .02), limits=c(0, 0.05)) +
      scale_fill_roadmap +
      facet_wrap(~ V1, ncol=1, scales='free') +
      theme_nature +
      theme(axis.text.x=element_blank(),
            panel.margin=unit(2, 'mm'),
            legend.position='none'))
Cairo(type="pdf", file=sub(".txt", ".pdf", args[1]), dpi='auto', width=183, height=120, unit="mm")
print(p)
dev.off()

warnings()
