library(graphite)
library(biomaRt)
library(RCurl)
library(stringr)
library(jsonlite)
library(httr)

#import filename
args=(commandArgs(TRUE))
if(length(args)==0){
  print("No arguments supplied.")
}
#setwd("~/Desktop")
infile <- args[1]
#infile <- "edgeR_dge_results_treatmentRSV_A549-treatmentmock_RSV_A549.txt_entrez.tsv_2020-04-07_06-11-06_SPIA_Results.csv"
outfile <- paste0(infile,"-Treatments.tsv")
setwd("~/fsl_groups/fslg_PickettLabGroup/spia")
#setwd("~/")

df_init <- data.frame(matrix(ncol = 9, nrow = 0))
colnames(df_init) <- c("Symbol","Name","Ensembl_ID","Modulation","Disorder","Molecule_Name","Molecule_Type","Database","Pathway_Name")
#y <- as.data.frame(c("1","2","3","4","5","6","7","8","9","10","11"))
write.table(df_init, file = outfile, row.names = FALSE, col.names=TRUE, sep = "\t", append = FALSE)

#import relevant pathways
humanReactome <- pathways("hsapiens", "reactome")
humanKEGG <- pathways("hsapiens", "kegg")
humanBioCarta <- pathways("hsapiens", "biocarta")
humanNCI <- pathways("hsapiens", "nci")
humanPanther <- pathways("hsapiens", "panther")

#read in spia results table
#setwd("~/Downloads")
in_data <- read.csv(file = infile, header = TRUE, sep = ",")
sig_paths <- as.vector(in_data$Name)
sig_dbs <- as.vector(in_data$SourceDB)

for(i in 1:length(sig_dbs)){
  #i <- 1
  print(paste0("Working...Pathway: ",sig_paths[i]))
  if(length(humanKEGG@entries[[sig_paths[i]]])==0){
    print(paste0("Skipping...Pathway: ",sig_paths[i]," -- No Protein Data"))
    next()
  }
  if(sig_dbs[i] == "KEGG"){
    #pg <- pathwayGraph(humanKEGG$'Cytokine-cytokine receptor interaction')
    db_entrezID <- as.vector(humanKEGG@entries[[sig_paths[i]]]@protEdges[["src"]])
    db_entrezID <- c(db_entrezID,humanKEGG@entries[[sig_paths[i]]]@protEdges[["dest"]])
    db_entrezID <- db_entrezID[!duplicated(db_entrezID)]
  }else if(sig_dbs[i] == "Reactome"){
    db_entrezID <- as.vector(humanReactome@entries[[sig_paths[i]]]@protEdges[["src"]])
    db_entrezID <- c(db_entrezID,humanReactome@entries[[sig_paths[i]]]@protEdges[["dest"]])
    db_entrezID <- db_entrezID[!duplicated(db_entrezID)]
  }else if(sig_dbs[i] == "NCI"){
    db_entrezID <- as.vector(humanNCI@entries[[sig_paths[i]]]@protEdges[["src"]])
    db_entrezID <- c(db_entrezID,humanNCI@entries[[sig_paths[i]]]@protEdges[["dest"]])
    db_entrezID <- db_entrezID[!duplicated(db_entrezID)]
  }else if(sig_dbs[i] == "Panther"){
    db_entrezID <- as.vector(humanPanther@entries[[sig_paths[i]]]@protEdges[["src"]])
    db_entrezID <- c(db_entrezID,humanPanther@entries[[sig_paths[i]]]@protEdges[["dest"]])
    db_entrezID <- db_entrezID[!duplicated(db_entrezID)]
  }else if(sig_dbs[i] == "Biocarta"){
    db_entrezID <- as.vector(humanBioCarta@entries[[sig_paths[i]]]@protEdges[["src"]])
    db_entrezID <- c(db_entrezID,humanBioCarta@entries[[sig_paths[i]]]@protEdges[["dest"]])
    db_entrezID <- db_entrezID[!duplicated(db_entrezID)]
  }else{
    print(paste0("Error: no valid pathway: ",sig_dbs[i]))
    next;
  }
  #entrez.genes <- db_entrezID
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <- getBM(
    filters="entrezgene_id",
    attributes=c("ensembl_gene_id", "entrezgene_id"),#, "hgnc_symbol", "uniprot_gn_id"),#, "go_id", "name_1006"),
    values=db_entrezID,
    mart=mart,
    unique = TRUE)
  
  #remove duplicate entrezIDs
  genes <- genes[!duplicated(genes$entrezgene_id), ]
  ensembl_vector <- as.vector(genes$ensembl_gene_id)
  
  url3 <- "https://api.opentargets.io/v3/platform/public/evidence/filter?target="
  datatype3 <- "&datatype=target.gene_info.symbol,target.target_class,evidence.target2drug.action_type,drug"#disease.efo_info.label,unique_association_fields.chembl_molecules,
  
  #submit to openTargets.org
  for(k in 1:length(ensembl_vector)){
    #k <- 1
    target <- ensembl_vector[k]
    #target <- "ENSG00000185436"
    call3 <- URLencode(paste(url3,target,sep = ""))
    target_data <- fromJSON(call3)
    if(length(target_data[["data"]][["drug"]] > 0)){
      print(paste0("drugs found: ",ensembl_vector[k]))
      drug_merged <- as.data.frame(do.call("cbind", list(target_data[["data"]][["target"]][["gene_info"]][["symbol"]],target_data[["data"]][["target"]][["gene_info"]][["name"]],target_data[["data"]][["target"]][["gene_info"]][["geneid"]],target_data[["data"]][["target"]][["activity"]],target_data[["data"]][["disease"]][["name"]],target_data[["data"]][["drug"]][["molecule_name"]],target_data[["data"]][["drug"]][["molecule_type"]],sig_dbs[i],sig_paths[i])))
      #colnames(drug_merged) <- c("Symbol","Name","Ensembl_ID","Modulation","Disorder","Molecule_Name","Molecule_Type")
      write.table(drug_merged, file = outfile, row.names = FALSE, col.names=FALSE, sep = "\t", append = TRUE)
    }else{
      print(paste0("no drugs found: ",ensembl_vector[k]))
    }
  }
  
  ##submit to EBI API
  ##for(j in 1:length(ensembl_vector)){
  #  j <- 1
  #  options(stringsAsFactors = FALSE)
  #  url <- "https://www.ebi.ac.uk/chembl/api/data/target/search.json?q="
  #  #gene_id <- ensembl_vector[j]
  #  gene_id <- "ENSG00000136244"
  #  print(paste0("Iteration ",j," of ",length(ensembl_vector)," Gene: ",gene_id,""))
  #  call <- paste0(url,gene_id)
  #  #call
  #  results <- GET(call)
  #  #results
  #  results_text <- content(results,"text")
  #  results_json <- fromJSON(results_text, flatten = TRUE)
  #  if(length(results_json[["targets"]][["target_chembl_id"]])>0){
  #    chembl_id <- results_json[["targets"]][["target_chembl_id"]]
  #  }else{
  #    print(paste0("No target found for: ",gene_id,"."))
  #    next;
  #  }
  #  
  #  #lookup target in ChEMBL:
  #  url1 <- "https://www.ebi.ac.uk/chembl/target_report_card/"
  #  call1 <- paste0(url1,chembl_id)
  #  results1 <- GET(call1)
  #  results_text1 <- content(results1,"text")
  #  results_json1 <- fromJSON(results_text1, flatten = TRUE)
  
  #}
}
print("Done")  

  
  
