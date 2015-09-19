library(ggplot2)
library(Cairo)

source("~/code/enr/plot/color_roadmap.R")
source("~/code/enr/plot/theme_nature.R")

args <- commandArgs(TRUE)
d <- read.table(args[1], header=0)
d$V1 <- factor(d$V1, levels=enh_cluster_ordering)
p <- (ggplot(d, aes(x=V1, y=V2, ymin=(V2 - V3), ymax=(V2 + V3), color=V1)) +
      geom_pointrange(size=I(.25)) +
      scale_x_discrete(name="Enhancer module", breaks=levels(d$V1)) +
      scale_y_continuous(name=expression(h[g]^2), expand=c(0, 0)) +
      scale_color_roadmap +
      theme_nature +
      theme(axis.text.x=element_blank(),
            axis.title.y=element_text(angle=0),
            legend.position='none'))
Cairo(type="pdf", file=sub(".txt", ".pdf", args[1]), dpi='auto', width=183, height=50, unit="mm")
print(p)
dev.off()
