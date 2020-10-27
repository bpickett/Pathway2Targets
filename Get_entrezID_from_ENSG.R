library(biomaRt)
if(!require(tools)) install.packages("tools", repos = "http://cran.us.r-project.org")
#library(tools)
#library(SingleCellExperiment)
#library(SummarizedExperiment)

#define input data
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}
infile <- args[1]

setwd("~/fsl_groups/fslg_PickettLabGroup/ARMOR")
#top1 <- readRDS("tximeta_se.rds")
top2 <- read.csv(infile, sep = "\t")
#top2$gene_id <- rownames(top2)
ensembl.genes <- as.vector(top2$gene_id)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),#, "go"),
  values=ensembl.genes,
  mart=mart)
#class(genes)
#class(top2)

#rf_results <- read.csv("TCGA_100-importance.tsv", sep = "\t")
#rf_results$ensembl_genes <- rownames(rf_results)
#combined_results <- merge(genes,top2,by.x = "ensembl_gene_id", by.y = "ensembl_genes")
temp_out <- merge(genes,top2,by.x = "ensembl_gene_id", by.y = "gene_id")#for ARMOR
#write.table(combined_results,file="TCGA-100-importance.tsv",sep="\t")
temp_out$entrezid <- NULL

outfile <-paste0(tools::file_path_sans_ext(infile),"_entrez.tsv")
#write.table(temp_out,file=paste0(infile,"_entrez.tsv"),sep="\t", col.names=NA)
write.table(temp_out, file=outfile, sep="\t", col.names=NA)
#write.csv(temp_out,file="DEG_list.tsv", sep="\t")
