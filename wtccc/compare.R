# Compare WTCCC1 imputed p-values with our LD-expanded p-values
# Author: Abhishek Sarkar <aksarkar@mit.edu>
library(Cairo)
library(grid)
library(hexbin)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=TRUE, sep=' ')
p <- hexbin(d$imputed, d$expanded, xlab='Imputed p-value',
            ylab='Representative p-value', xbins=100)
CairoPNG(file='heatmap.png', height=1024, width=1024)
grid.newpage()
plot(p)
dev.off()
