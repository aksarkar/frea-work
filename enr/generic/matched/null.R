library(ggplot2)
d = read.table('trace')
p <- (qplot(x=d$V1, geom='density', xlab="Count of top resampled SNPs", ylab="Probability density") +
      geom_vline(xintercept=1721, color='red') +
      theme_bw())
pdf(height=4,width=6); print(p); dev.off()
