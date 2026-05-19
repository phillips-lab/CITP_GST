
library(edgeR)
library(openxlsx)
library(stringr)

# Adult ages to test.
EXP <- c(2,4,8,12,16)

dir.create("plots/mds", recursive=TRUE, showWarnings=FALSE)
dir.create("reports", recursive=TRUE, showWarnings=FALSE)

# Gene annotations.
genes <- readRDS("r_data/gene_annotations_Cel.Rdata")

# Load the feature counts.
fc.gst <- read.table("GST_featureCounts.txt", header=TRUE)
rownames(fc.gst) <- fc.gst$Geneid
fc.gst <- as.matrix(fc.gst[c(1:38)])

stopifnot(all(row.names(fc.gst)==genes$gene_id))

wb <- createWorkbook()

for(i in 1:length(EXP)) {
  
  fc <- fc.gst[,colnames(fc.gst) %in% colnames(fc.gst)[grepl(paste0("_",EXP[i],"_"),colnames(fc.gst))]]
  
  fc.col <- colnames(fc)
  col.n <- c()
  grp <- c()
  for(j in 1:length(fc.col)) {
    col.s <- str_split(fc.col[j], "_")[[1]]
    if(col.s[3]=="GST") { grp <- append(grp,"T") } else { grp <- append(grp,"C") }
    col.n <- append(col.n, paste0(col.s[3],"_d",col.s[2],"r",col.s[4]))
  }
  
  colnames(fc) <- col.n
  
  show(colnames(fc))
  
  # Grouping factor.
  group <- factor(grp)
  
  # Create the data object.
  y <- DGEList(counts=fc, group=group, genes=genes$gene_name)
  
  # Filter out lowly expressed genes.
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  
  # Normalize the library sizes.
  y <- calcNormFactors(y)
  
  # Save MDS plot.
  pdf(file=paste0("plots/mds/mds_day",EXP[i],".pdf"), width=8, height=8)
  plotMDS(y)
  dev.off()
  
  # Create design matrix.
  design <- model.matrix(~group)
  
  show(design)
  
  # Estimate common dispersion and tagwise dispersions.
  y <- estimateDisp(y, design)
  
  # Perform quasi-likelihood F-tests for DE.
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit)
  
  # Save table.
  tags <- as.data.frame(topTags(qlf, n=nrow(qlf)))
  tags <- cbind(row.names(tags), tags)
  colnames(tags) <- c("geneWbid","gene","logFC","logCPM","F","pValue","FDR")
  
  pf <- c()
  for(j in 1:nrow(tags)) {
    pf <- append(pf, ifelse(tags$FDR[j]<=.05 && abs(tags$logFC[j])>=1, TRUE, FALSE))
  }
  tags$passFilter <- pf
  
  addWorksheet(wb=wb, sheetName=paste0("de_day",EXP[i]))
  modifyBaseFont(wb, fontSize=16)
  writeDataTable(wb=wb, sheet=i, x=tags)
  
}

saveWorkbook(wb, paste0("reports/two_group_comparisons.xlsx"), overwrite=TRUE)
