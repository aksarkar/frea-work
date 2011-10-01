args <- commandArgs(TRUE)
df <- read.delim(args[1], header=0)
print(wilcox.test(subset(df, V6 == 1)$V5, subset(df, V6 == 0)$V5, alternative='less')$p.value)
