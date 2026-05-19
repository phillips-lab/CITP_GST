
library(plyranges)

GTF_FN <- "gtf/Caenorhabditis_elegans.WBcel235.112.gtf"

dir.create("r_data", recursive=TRUE, showWarnings=FALSE)

gff <- read_gff(GTF_FN) %>% select(gene_id, gene_name)

genes <- unique(data.frame(gene_id=gff$gene_id, gene_name=gff$gene_name))

saveRDS(genes, file="r_data/gene_annotations_Cel.Rdata")
