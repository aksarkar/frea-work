library(directlabels)
library(ggplot2)
library(Cairo)

args <- commandArgs(TRUE)
d <- read.table(args[1], header=0)
m <- lm(V4 ~ V2, data=d)
print(summary(m))
e <- fortify(m, d)
e$top <- e$V1 %in% c("CD4+_CD25int_CD127+_Tmem_Primary_Cells",
"CD4+_CD25-_CD45RO+_Memory_Primary_Cells",
"CD56_Primary_Cells",
"CD3_Primary_Cells_Peripheral_UW",
"CD4+_CD25-_Th_Primary_Cells",
"CD8_Memory_Primary_Cells",
"CD8_Naive_Primary_Cells",
"CD4+_CD25-_CD45RA+_Naive_Primary_Cells",
"CD4_Naive_Primary_Cells",
"CD4+_CD25+_CD127-_Treg_Primary_Cells") ## ,
## "CD4_Memory_Primary_Cells",
## "CD3_Primary_Cells_Cord_BI",
## "Peripheral_Blood_Mononuclear_Primary_Cells",
## "CD4+_CD25-_IL17+_PMA-Ionomcyin_stimulated_Th17_Primary_Cells",
## "CD4+_CD25-_IL17-_PMA-Ionomycin_stimulated_MACS_purified_Th_Primary_Cells")
p <- (ggplot(e, aes(x=V2, y=V4, color=top)) +
      scale_x_continuous(name="Number of SNPs tagging top T1D—associated enhancer") +
      scale_y_continuous(name=expression(h[g]^2), limits=c(0, .25)) +
      scale_color_manual(values=c('TRUE' = 'red', 'FALSE' = 'black'), guide='none') +
      stat_smooth(method='lm', se=0, color='red') +
      ## geom_dl(data=F, aes(label=V1), method='smart.grid') +
      geom_point() +
      theme_bw() +
      theme(axis.text.x=element_text(angle=-90, hjust=0, vjust=.5),
            axis.title.y=element_text(angle=0),
            panel.grid.minor=element_blank()))
Cairo(type="svg", name=sub(".in", ".svg", args[1]), dpi=96, width=267, height=167, unit="mm")
print(p)
dev.off()
