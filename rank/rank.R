library(Cairo)
library(ggplot2)
source("/home/unix/aksarkar/code/enr/plot/theme_nature.R")
args <- commandArgs(TRUE)
d <- read.table(args[1])
p <- (qplot(data=d, x=V1, xlab='Top n SNPs (full meta-analysis)', y=V3, ylab='Hold-out Pearson correlation',
            color=factor(V2), geom='line', size=I(.25)) +
      geom_hline(yintercept=0, color="black", size=I(.1), shape='dashed') +
      scale_x_continuous(limits=c(0, 1e5)) +
      scale_color_brewer(palette="Dark2", name="Cohort") +
      theme_nature +
      theme(legend.position="right"))
Cairo(file=sub('.txt$', '.pdf', args[1]), type='pdf', height=50, width=89, unit='mm')
print(p)
dev.off()
