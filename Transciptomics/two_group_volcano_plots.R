
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(openxlsx)

# Adult ages to plot.
EXP <- c(8,12,16)

dir.create("plots/volcano", recursive=TRUE, showWarnings=FALSE)

for(i in 1:length(EXP)) {
  
  # Load the DE results for the adult age.
  de <- read.xlsx("reports/two_group_comparisons.xlsx", sheet=i)
  
  de$delabel <- NA
  de$delabel[de$FDR<=.05 & abs(de$logFC)>=1] <- de$gene[de$FDR<=.05 & abs(de$logFC)>=1]
  de$diffTrend <- "no"
  de$diffTrend[de$FDR<=.05 & de$logFC>=1] <- "up"
  de$diffTrend[de$FDR<=.05 & de$logFC<=-1] <- "down"
  
  p <- ggplot(data=de, aes(x=logFC, y=-log10(FDR), col=diffTrend, label=delabel)) +
     geom_vline(xintercept=c(-1,1), linetype="dotted", linewidth=0.8, col="#222222") +
     geom_hline(yintercept=-log10(0.05), linetype="dotted", linewidth=0.8, col="#222222") +
     geom_point() +
     theme_classic() +
     theme(legend.position="none", plot.title=element_text(size=18, face="bold", hjust=0.5)) +
     theme(axis.line=element_line(size=1), axis.text=element_text(size=16)) +
     geom_text_repel(size=5) +
     annotate(geom="label", x=-5, y=7, label="Downregulated", size=5.2, color="#2F67B1") +
     annotate(geom="label", x=9, y=7, label="Upregulated", size=5.2, color="#BF2C23") +
     scale_color_manual(values=c("#2F67B1","#444444","#BF2C23")) +
     scale_x_continuous(limits=c(-6,10), breaks=c(-5,-3,-1,1,3,5,7,9)) +
     scale_y_continuous(limits=c(0,7), breaks=c(1,3,5,7)) +
     ggtitle(paste("Day",EXP[i])) +
     labs(x=expression(log[2]~Fold~Change), y=expression(-log[1][0]~FDR))
  
  ggsave(paste0("plots/volcano/volcano_plot_day",EXP[i],".png"), width=10, height=5, p, device="png")
  
}
