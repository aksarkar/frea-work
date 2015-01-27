library(reshape)

args <- commandArgs(TRUE)
data <- read.table(args[1], sep=' ', header=FALSE)
samples <- read.delim('/broad/compbio/aksarkar/data/roadmap/meta/sample_info_WM20140407.txt')
data$V3 <- factor(data$V3, levels=samples$EID[samples$position])
data <- transform(data, z=(V5 - V6) / sqrt(V7))
matrix <- cast(data, V3 ~ V1, value='z', fun.aggregate=max)
write.table(matrix, file=sub('match', 'permutation-test', args[1]), quote=FALSE, sep='\t', row.names=FALSE)
