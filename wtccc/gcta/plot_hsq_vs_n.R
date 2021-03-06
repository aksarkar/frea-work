library(ggplot2)
library(grid)
library(plyr)
library(Cairo)

source("~/code/enr/plot/theme_nature.R")

jacknife_se <- function(x, label) {
    x$w <- 1 / x$V4 ** 2
    data.frame(mean=mean(x$V3), se=sd(x$V3), order=rep(label))
}

args <- commandArgs(TRUE)
hsq_by_n <- ddply(read.table(args[1], header=0, sep=' '), .(V2, V5), jacknife_se, "Observed")
hsq_by_random_n <- ddply(read.table(args[2], header=0, sep=' '), .(V2, V5), jacknife_se, "Random")
hsq <- rbind(hsq_by_n, hsq_by_random_n)
## all_sample_hsq <- read.table(args[3], header=0, sep=' ')
p <- (ggplot(hsq, aes(x=V2, y=mean, ymin=mean - se, ymax=mean + se)) +
      geom_line(size=I(.25), aes(color=order)) +
      geom_ribbon(alpha=I(.2), aes(fill=order)) +
      ## geom_hline(data=all_sample_hsq, aes(yintercept=V3), size=I(.25), linetype='dashed') +
      geom_hline(yintercept=0, size=I(.25)) +
      scale_x_log10(name='Top n tags', labels=unique(hsq$V2), breaks=unique(hsq$V2)) +
      scale_y_continuous(name=expression(h[g]^2)) +
      scale_color_manual(values=c("Observed"="red", "Random"="gray50"), name="Ordering") +
      scale_fill_manual(values=c("Observed"="red", "Random"="gray50"), name="Ordering") +
      facet_wrap(~ V5, ncol=4, scales='free') +
      theme_nature +
      theme(axis.title.y=element_text(angle=0),
            axis.text=element_text(size=5),
            axis.text.x=element_text(angle=-90, hjust=0, vjust=.5),
            legend.position="right",
            panel.margin=unit(rep(1, 4), 'mm')
            ))
Cairo(dpi='auto', file=sub('.txt', '.pdf', args[1]), height=80, type='pdf', units='mm', width=210)
print(p)
dev.off()
