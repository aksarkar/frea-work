library(Cairo)
library(directlabels)
library(ggplot2)
library(gridExtra)
library(plyr)
library(scales)

source("~/code/enr/plot/color_roadmap.R")
source("~/code/enr/plot/theme_nature.R")

heatmap <- function(data) {
  ggplotGrob(ggplot(data, aes(x=V1, y=V3, fill=fold)) +
             geom_tile() +
             scale_x_discrete(limits=c('iibdgc-cd', 'stahl-ra', 'wtccc-t1d',
                                'diagram-t2d', 'cardiogram-cad', 'igap-ad',
                                'pgc-bip', 'pgc-scz'),
                              labels=c('CD', 'RA', 'T1D', 'T2D', 'CAD', 'AD',
                                'BD', 'SZ')) +
             coord_fixed() +
             facet_grid(. ~ group_name) +
             heatmap_ramp() +
             theme_nature +
             theme(text=element_text(size=5),
                   axis.title=element_blank(),
                   axis.text.x=element_text(angle=-90, hjust=0, vjust=.5),
                   strip.text=element_text(size=5),
                   strip.text.y=element_text(angle=0, hjust=0)))
}


args <- commandArgs(TRUE)
data <- read.table(args[1], header=FALSE, sep=' ')
data <- transform(data, fold=V5 / V6)
data <- subset(data, V8 < .05 / nrow(data))
samples <- read.delim('/broad/compbio/aksarkar/data/roadmap/meta/sample_info_WM20140407.txt')
data <- merge(x=data, y=samples, by.x='V3', by.y='EID')
data$V3 <- factor(data$V3, levels=rev(samples$EID[samples$position]))
ps <- dlply(data, .(group_name), heatmap)
grobs <- do.call(rbind, c(ps, size='last'))
Cairo(type='pdf', dpi='auto', width=189, height=400, units='mm')
grid.draw(grobs)
dev.off()
