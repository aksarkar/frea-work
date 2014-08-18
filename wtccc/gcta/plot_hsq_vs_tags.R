library(directlabels)
library(ggplot2)
library(Cairo)

source("~/code/enr/plot/color_roadmap.R")

args <- commandArgs(TRUE)
d <- read.table(args[1], header=0)
m <- lm(V4 ~ V2, data=d)
e <- fortify(m, d)
if (length(args) > 1) {
    p <- ggplot(e, aes(x=V2, y=V4, color=factor(V1))) + scale_color_roadmap
} else {
    p <- ggplot(e, aes(x=V2, y=V4))
}
p <- (p +
      scale_x_continuous(name="Number of SNPs in top associated enhancers") +
      scale_y_continuous(name=expression(h[g]^2)) +
      stat_smooth(method='lm', se=0, color='red') +
      geom_point() +
      theme_bw() +
      theme(axis.text.x=element_text(angle=-90, hjust=0, vjust=.5),
            axis.title.y=element_text(angle=0),
            panel.grid.minor=element_blank(),
            legend.position='none'))
Cairo(type="pdf", name=sub(".in", ".pdf", args[1]), dpi=96, width=160, height=100, unit="mm")
print(p)
dev.off()
