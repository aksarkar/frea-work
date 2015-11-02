library(ggplot2)
library(grid)
library(plyr)
library(Cairo)

source("~/code/enr/plot/theme_nature.R")

inv_var <- function(x, label) {
    x$w <- 1 / x$V4 ** 2
    data.frame(mean=sum(x$V3 * x$w) / sum(x$w), se=sqrt(1 / sum(x$w)), order=rep(label))
}

args <- commandArgs(TRUE)
hsq_by_n <- ddply(read.table(args[1], header=0, sep=' '), .(V1, V5), inv_var, "Observed")
hsq_by_random_n <- ddply(read.table(args[2], header=0, sep=' '), .(V1, V5), inv_var, "Random")
hsq <- rbind(hsq_by_n, hsq_by_random_n)
all_sample_hsq <- read.table(args[3], header=0, sep=' ')
p <- (ggplot(hsq, aes(x=V1, y=mean, ymin=mean - se, ymax=mean + se, color=order)) +
      geom_line(size=I(.25)) +
      geom_ribbon(alpha=I(.25)) +
      geom_hline(data=all_sample_hsq, aes(yintercept=V3), size=I(.25), linetype='dashed') +
      scale_x_log10(name='Top n tags', labels=unique(hsq$V1), breaks=unique(hsq$V1)) +
      scale_y_continuous(name=expression(h[g]^2)) +
      scale_color_manual(values=c("Observed"="black", "Random"="red")) +
      facet_wrap(~ V5, ncol=2, scales='free') +
      theme_nature +
      theme(axis.title.y=element_text(angle=0),
            axis.text=element_text(size=5),
            axis.text.x=element_text(angle=-90, hjust=0, vjust=.5),
            legend.position="right",
            panel.margin=unit(rep(5, 4), 'mm'),
            plot.margin=unit(c(1, 0, 0, 0), 'mm')
            ))
Cairo(dpi='auto', file=sub('.txt', '.pdf', args[1]), height=190, type='pdf', units='mm', width=160)
print(p)
dev.off()
