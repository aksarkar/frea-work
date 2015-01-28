library(ggplot2)
library(Cairo)

source("~/code/enr/plot/theme_nature.R")

marchenko.pastur <- function(x, lambda) {
  upper <- (1 + sqrt(lambda)) ** 2
  lower <- (1 - sqrt(lambda)) ** 2
  sqrt(as.numeric(x >= lower) * as.numeric(x <= upper) * (upper - x) * (x - lower)) / (2 * pi * lambda * x)
}

args <- commandArgs(TRUE)
eigenvals <- read.delim(args[1], header=FALSE)
lambda <- as.numeric(args[2])
Cairo(type='pdf', width=83, height=83, units='mm')
(ggplot(eigenvals, aes(x=V1)) +
 geom_density(size=I(.35 / ggplot2:::.pt)) +
 stat_function(fun=marchenko.pastur, args=list(lambda=lambda), geom='path', size=I(.35 / ggplot2:::.pt), color='red') +
 theme_nature)
dev.off()

print(eigenvals[eigenvals$V1 > (1 + sqrt(lambda)) ** 2,])
