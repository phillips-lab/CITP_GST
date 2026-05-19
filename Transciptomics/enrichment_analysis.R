
library(clusterProfiler)
library(openxlsx)
library(org.Ce.eg.db)
library(ReactomePA)
library(tidyverse)

de <- list()

# Load DE restuls for adult ages 8, 12, 16.
de[[1]] <- read.xlsx("reports/two_group_comparisons.xlsx", sheet=3)
de[[2]] <- read.xlsx("reports/two_group_comparisons.xlsx", sheet=4)
de[[3]] <- read.xlsx("reports/two_group_comparisons.xlsx", sheet=5)

ages <- c("Day 8","Day 12","Day 16")

directions <- c("both","up","down")

# Convert from WormBase to Entrez IDs.
wb2eg <- function(genes) {
  bitr(genes, fromType="WORMBASE", toType="ENTREZID", OrgDb="org.Ce.eg.db")
}

# Gene universe.
genes <- unique(c(de[[1]]$geneWbid, de[[2]]$geneWbid, de[[3]]$geneWbid))
uni.eg <- wb2eg(genes)

gobp.all <- data.frame()
gomf.all <- data.frame()
reac.all <- data.frame()
kegg.all <- data.frame()

# Perform enrichments at each age.
for(i in 1:length(ages)) {
  
  cat(paste("Processing",ages[i],"\n"))
  
  gobp.df <- data.frame()
  gomf.df <- data.frame()
  reac.df <- data.frame()
  kegg.df <- data.frame()
  
  eg <- list()
  
  # all DE genes
  g <- subset(de[[i]], abs(logFC)>=1 & FDR<=.05)
  if(any(g$geneWbid %in% uni.eg$WORMBASE)) {
    eg[[1]] <- wb2eg(g$geneWbid)
  } else {
    eg[[1]] <- data.frame()
  }
  
  # upregulated genes
  g <- subset(de[[i]], logFC>=1 & FDR<=.05)
  if(any(g$geneWbid %in% uni.eg$WORMBASE)) {
    eg[[2]] <- wb2eg(g$geneWbid)
  } else {
    eg[[2]] <- data.frame()
  }
  
  # downregulated genes
  g <- subset(de[[i]], logFC<=-1 & FDR<=.05)
  if(any(g$geneWbid %in% uni.eg$WORMBASE)) {
    eg[[3]] <- wb2eg(g$geneWbid)
  } else {
    eg[[3]] <- data.frame()
  }
  
  for(j in 1:length(directions)) {
    
    if(nrow(eg[[j]]) > 0) {
      
      # GOBP enrichment.
      en.gobp <- enrichGO(gene=eg[[j]]$ENTREZID,
                            OrgDb=org.Ce.eg.db,
                            ont="BP",
                            universe=uni.eg$ENTREZID,
                            readable=TRUE)
      
      en.gobp <- as.data.frame(en.gobp)
      if(nrow(en.gobp)>0) {
        en.gobp$Age <- ages[i]
        en.gobp$Direction <- directions[j]
      }
      gobp.df <- rbind(gobp.df, en.gobp)
      
      # GOMF enrichment.
      en.gomf <- enrichGO(gene=eg[[j]]$ENTREZID,
                          OrgDb=org.Ce.eg.db,
                          ont="MF",
                          universe=uni.eg$ENTREZID,
                          readable=TRUE)
      
      en.gomf <- as.data.frame(en.gomf)
      if(nrow(en.gomf)>0) {
        en.gomf$Age <- ages[i]
        en.gomf$Direction <- directions[j]
      }
      gomf.df <- rbind(gomf.df, en.gomf)
      
      # Reactome enrichment.
      en.reac <- enrichPathway(eg[[j]]$ENTREZID,
                               organism="celegans",
                               universe=uni.eg$ENTREZID)
      
      en.reac <- as.data.frame(en.reac)
      if(nrow(en.reac)>0) {
        en.reac$Age <- ages[i]
        en.reac$Direction <- directions[j]
      }
      reac.df <- rbind(reac.df, en.reac)
      
      # KEGG enrichment.
      en.kegg <- enrichKEGG(gene=eg[[j]]$ENTREZID,
                            keyType="ncbi-geneid",
                            organism="cel",
                            universe=uni.eg$ENTREZID)
      
      en.kegg <- as.data.frame(en.kegg)
      if(nrow(en.kegg)>0) {
        en.kegg$Age <- ages[i]
        en.kegg$Direction <- directions[j]
      }
      kegg.df <- rbind(kegg.df, en.kegg)
      
    }
    
  }
  
  if(nrow(gobp.df)>0) { gobp.all <- rbind(gobp.all, gobp.df[c(10:11,1:9)]) }
  if(nrow(gomf.df)>0) { gomf.all <- rbind(gomf.all, gomf.df[c(10:11,1:9)]) }
  if(nrow(reac.df)>0) { reac.all <- rbind(reac.all, reac.df[c(10:11,1:9)]) }
  if(nrow(kegg.df)>0) { kegg.all <- rbind(kegg.all, kegg.df[c(12:13,1:11)]) }
  
}

wb <- createWorkbook()
addWorksheet(wb=wb, sheetName="GO Biological Process")
writeDataTable(wb=wb, sheet=1, x=gobp.all)
addWorksheet(wb=wb, sheetName="GO Molecular Function")
writeDataTable(wb=wb, sheet=2, x=gomf.all)
addWorksheet(wb=wb, sheetName="Reactome Pathway")
writeDataTable(wb=wb, sheet=3, x=reac.all)
addWorksheet(wb=wb, sheetName="KEGG Pathway")
writeDataTable(wb=wb, sheet=4, x=kegg.all)
saveWorkbook(wb, "reports/enrichment_analysis.xlsx", overwrite=TRUE)
