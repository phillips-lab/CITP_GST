
library(ggplot2)
library(openxlsx)
library(tidyverse)

dir.create("plots/enrichments", recursive=TRUE, showWarnings=FALSE)

AGES <- c("Day 8","Day 12","Day 16")
ENS <- c("GO Biological Processes","GO Molecular Functions","Reactome Pathways","KEGG Pathways")
ENS_ABBR <- c("gobp","gomf","reactome","kegg")
LABELS <- c("8","12","16")

make.plot <- function(dat, enrich, width, height, dir, dir.text) {
  
  gene.ratio <- c()
  for(j in 1:nrow(dat)) {
    gene.ratio <- append(gene.ratio, as.numeric(strsplit(dat$GeneRatio[j],"/")[[1]][1])/as.numeric(strsplit(dat$GeneRatio[j],"/")[[1]][2]))
  }
  
  dat$`Gene ratio` <- gene.ratio
  dat$Significance <- -log10(dat$p.adjust)
  
  # insert line breaks in long term names
  dat$Description[dat$Description=="biological process involved in interspecies interaction between organisms"] <- "biological process involved in interspecies\ninteraction between organisms"
  dat$Description[dat$Description=="oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen"] <- "oxidoreductase activity, acting on paired donors,\nwith incorporation or reduction of molecular oxygen"
  dat$Description[dat$Description=="oxidoreductase activity, acting on paired donors, with incorporation or reduction of molecular oxygen, reduced flavin or flavoprotein as one donor, and incorporation of one atom of oxygen"] <- "oxidoreductase activity, acting on paired donors,\nwith incorporation or reduction of molecular oxygen,\nreduced flavin or flavoprotein as one donor,\nand incorporation of one atom of oxygen"
  dat$Description[dat$Description=="RNA polymerase II transcription regulatory region sequence-specific DNA binding"] <- "RNA polymerase II transcription regulatory\nregion sequence-specific DNA binding"
  dat$Description[dat$Description=="positive regulation of nucleobase-containing compound metabolic process"] <- "positive regulation of nucleobase-containing\ncompound metabolic process"
  dat$Description[dat$Description=="SCF-dependent proteasomal ubiquitin-dependent protein catabolic process"] <- "SCF-dependent proteasomal ubiquitin-dependent\nprotein catabolic process"
  dat$Description[dat$Description=="DNA-binding transcription factor activity, RNA polymerase II-specific"] <- "DNA-binding transcription factor activity,\nRNA polymerase II-specific"
  dat$Description[dat$Description=="RNA polymerase II cis-regulatory region sequence-specific DNA binding"] <- "RNA polymerase II cis-regulatory region\nsequence-specific DNA binding"
  dat$Description[dat$Description=="RUNX1 regulates genes involved in megakaryocyte differentiation and platelet function"] <- "RUNX1 regulates genes involved in megakaryocyte\ndifferentiation and platelet function"
  dat$Description[dat$Description=="Synthesis of (16-20)-hydroxyeicosatetraenoic acids (HETE)"] <- "Synthesis of (16-20)-\nhydroxyeicosatetraenoic acids (HETE)"
  
  p <- ggplot(data=dat,
              aes(x=factor(Age, level=AGES), y=forcats::fct_rev(factor(Description)), color=Significance, size=`Gene ratio`)) +
       geom_point() +
       scale_color_gradient(low="blue", high="red") +
       xlab("Age in days of adulthood") +
       ylab("") +
       ggtitle(paste0(ENS[i], dir.text)) +
       scale_x_discrete(labels=LABELS) +
       theme_bw() +
       theme(plot.title=element_text(size=12), legend.title=element_text(size=9), panel.border=element_rect(colour="#333333", fill=NA, linewidth=.7))
  
  ggsave(paste0("plots/enrichments/enrichment_plot_",dir,"_",enrich,".png"), width=width, height=height, p, device="png")
  
}


# All DE genes.

WIDTH  <- c(4.5, 5, 5, 4.7)
HEIGHT <- c(4, 5.5, 5.1, 4.8)

for(i in 1:length(ENS)) {
  
  en <- read.xlsx("reports/enrichment_analysis.xlsx", sheet=i)
  en <- subset(en, Age %in% AGES & Direction=="both")
  dat <- data.frame()
  
  # obtain the top ten hits at each age
  for(j in 1:length(AGES)) {
    df <- subset(en, Age==AGES[j])
    if(nrow(df)<10) { top<-nrow(df) } else { top<-10 }
    dat <- rbind(dat, df[1:top,])
  }
  
  terms <- unique(dat$Description)
  dat <- subset(en, Description %in% terms)
  
  p <- make.plot(dat, ENS_ABBR[i], WIDTH[i], HEIGHT[i], "all_deg", "\nAll DE genes")
  
}


# Upregulated genes.

WIDTH  <- c(4.6, 5.1, 4.2, 4.9)
HEIGHT <- c(4, 4.8, 4.7, 4.7)

for(i in 1:length(ENS)) {
  
  en <- read.xlsx("reports/enrichment_analysis.xlsx", sheet=i)
  en <- subset(en, Age %in% AGES & Direction=="up")
  dat <- data.frame()
  
  # obtain the top ten hits at each age
  for(j in 1:length(AGES)) {
    df <- subset(en, Age==AGES[j])
    if(nrow(df)<10) { top<-nrow(df) } else { top<-10 }
    dat <- rbind(dat, df[1:top,])
  }
  
  terms <- unique(dat$Description)
  dat <- subset(en, Description %in% terms)
  
  p <- make.plot(dat, ENS_ABBR[i], WIDTH[i], HEIGHT[i], "upreg", "\nUpregulated genes")
  
}


# Downregulated genes.

WIDTH  <- c(5.6, 5.2, 5.1, 4.6)
HEIGHT <- c(6.1, 5, 4.1, 3.5)

for(i in 1:length(ENS)) {
  
  en <- read.xlsx("reports/enrichment_analysis.xlsx", sheet=i)
  en <- subset(en, Age %in% AGES & Direction=="down")
  dat <- data.frame()
  
  # obtain the top ten hits at each age
  for(j in 1:length(AGES)) {
    df <- subset(en, Age==AGES[j])
    if(nrow(df)<10) { top<-nrow(df) } else { top<-10 }
    dat <- rbind(dat, df[1:top,])
  }
  
  terms <- unique(dat$Description)
  dat <- subset(en, Description %in% terms)
  
  p <- make.plot(dat, ENS_ABBR[i], WIDTH[i], HEIGHT[i], "downreg", "\nDownregulated genes")
  
}
