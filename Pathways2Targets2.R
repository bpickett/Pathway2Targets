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
#infile <- "edgeR_dge_results_treatmentRSV_A549-treatmentmock_RSV_A549.txt_entrez.tsv_2020-05-03_12-57-17_SPIA_Results.csv"
outfile <- paste0(infile,"-Treatments.tsv")
setwd("~/fsl_groups/fslg_PickettLabGroup/spia")
#setwd("~/")
merged_drugs <- as.data.frame(NULL)

#df_init <- data.frame(matrix(ncol = 11, nrow = 0))
#colnames(df_init) <- c("Target_ID","Target_Symbol","Target_Name","Drug_ID","Drug_Name","Is_FDA_Approved","Highest_Clinical_Trial_Phase","Has_Been_Withdrawn","Approved_Indications","Pathway_DB","Pathway_Name")

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
  print(paste0("Working...Pathway ",i," of ", length(sig_dbs)," : ",sig_paths[i]))
  
  if(sig_dbs[i] == "KEGG"){
    if(length(humanKEGG@entries[[sig_paths[i]]])==0){
     print(paste0("Skipping...Pathway ",i,"of ", length(sig_dbs)," : ",sig_paths[i]," -- No Protein Data"))
     next()
    }
    db_entrezID <- as.vector(humanKEGG@entries[[sig_paths[i]]]@protEdges[["src"]])
    db_entrezID <- c(db_entrezID,humanKEGG@entries[[sig_paths[i]]]@protEdges[["dest"]])
    db_entrezID <- db_entrezID[!duplicated(db_entrezID)]
    
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    genes <- getBM(
      filters="entrezgene_id",
      attributes=c("ensembl_gene_id", "entrezgene_id"),#, "hgnc_symbol", "uniprot_gn_id"),#, "go_id", "name_1006"),
      values=db_entrezID,
      mart=mart,
      unique = TRUE)
  }else if(sig_dbs[i] == "Reactome"){
    if(length(humanReactome@entries[[sig_paths[i]]])==0){
      print(paste0("Skipping...Pathway ",i,"of ", length(sig_dbs)," : ",sig_paths[i]," -- No Protein Data"))
      next()
    }
    db_entrezID <- as.vector(humanReactome@entries[[sig_paths[i]]]@protEdges[["src"]])
    db_entrezID <- c(db_entrezID,humanReactome@entries[[sig_paths[i]]]@protEdges[["dest"]])
    db_entrezID <- db_entrezID[!duplicated(db_entrezID)]
    
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    genes <- getBM(
      filters="uniprotswissprot",
      attributes=c("uniprotswissprot", "ensembl_gene_id", "entrezgene_id"),#, "hgnc_symbol", "uniprot_gn_id"),#, "go_id", "name_1006"),
      values=db_entrezID,
      mart=mart,
      unique = TRUE)
  }else if(sig_dbs[i] == "NCI"){
    if(length(humanNCI@entries[[sig_paths[i]]])==0){
      print(paste0("Skipping...Pathway ",i,"of ", length(sig_dbs)," : ",sig_paths[i]," -- No Protein Data"))
      next()
    }
    db_entrezID <- as.vector(humanNCI@entries[[sig_paths[i]]]@protEdges[["src"]])
    db_entrezID <- c(db_entrezID,humanNCI@entries[[sig_paths[i]]]@protEdges[["dest"]])
    db_entrezID <- db_entrezID[!duplicated(db_entrezID)]
    
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    genes <- getBM(
      filters="uniprotswissprot",
      attributes=c("uniprotswissprot", "ensembl_gene_id", "entrezgene_id"),#, "hgnc_symbol", "uniprot_gn_id"),#, "go_id", "name_1006"),
      values=db_entrezID,
      mart=mart,
      unique = TRUE)
  }else if(sig_dbs[i] == "Panther"){
    if(length(humanPanther@entries[[sig_paths[i]]])==0){
      print(paste0("Skipping...Pathway ",i,"of ", length(sig_dbs)," : ",sig_paths[i]," -- No Protein Data"))
      next()
    }
    db_entrezID <- as.vector(humanPanther@entries[[sig_paths[i]]]@protEdges[["src"]])
    db_entrezID <- c(db_entrezID,humanPanther@entries[[sig_paths[i]]]@protEdges[["dest"]])
    db_entrezID <- db_entrezID[!duplicated(db_entrezID)]
    
    mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
    genes <- getBM(
      filters="uniprotswissprot",
      attributes=c("uniprotswissprot", "ensembl_gene_id", "entrezgene_id"),#, "hgnc_symbol", "uniprot_gn_id"),#, "go_id", "name_1006"),
      values=db_entrezID,
      mart=mart,
      unique = TRUE)
  }else if(sig_dbs[i] == "Biocarta"){
    if(length(humanBioCarta@entries[[sig_paths[i]]])==0){
      print(paste0("Skipping...Pathway ",i,"of ", length(sig_dbs)," : ",sig_paths[i]," -- No Protein Data"))
      next()
    }
    db_entrezID <- as.vector(humanBioCarta@entries[[sig_paths[i]]]@protEdges[["src"]])
    db_entrezID <- c(db_entrezID,humanBioCarta@entries[[sig_paths[i]]]@protEdges[["dest"]])
    db_entrezID <- db_entrezID[!duplicated(db_entrezID)]
  }else{
    print(paste0("Error: no valid pathway: ",sig_dbs[i]))
    next;
  }
  entrez.genes <- db_entrezID
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
  
  #build query string:
  #**note, does not use "tradeNames" or "approvedIndications" as data under "drug" table
  query_string = "
       query target($ensemblId: String!){
    target(ensemblId: $ensemblId){
      id
      approvedSymbol
      approvedName
    knownDrugs {
      uniqueDrugs
      rows {
        drug {
          id
          name
          isApproved
          maximumClinicalTrialPhase
          hasBeenWithdrawn
        }
      }
    }
   }
  }
  "
  # Set base URL of GraphQL API endpoint
  base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
  
  #submit to openTargets.org
  for(k in 1:length(ensembl_vector)){
    #k <- 5
    target <- ensembl_vector[k]
    #target <- "ENSG00000232810"
    #target <- "ENSG00000197919"
    #target <- "ENSG00000185436"
    # Set variables object of arguments to be passed to endpoint
    variables <- list("ensemblId" = target)
    # Construct POST request body object with query string and variables
    post_body <- list(query = query_string, variables = variables)
    # Perform POST request
    r <- POST(url=base_url, body=post_body, encode='json')
    # print to console
    #print(content(r)$data)
    target_data <- content(r)$data
    
    #if(length(target_data[["data"]][["drug"]] > 0)){
    if(length(unlist(target_data[["target"]][["knownDrugs"]])) == 0){
      print(paste0("Pathway Member ",k," of ",length(ensembl_vector)," No drugs found: ",ensembl_vector[k]))
    }else{
      print(paste0("Pathway Member ",k," of ",length(ensembl_vector)," Drug(s) found:  ", ensembl_vector[k]))
      
      for(y in 1:length(target_data[["target"]][["knownDrugs"]][["rows"]])){
        #y <- 1
        unique_chembl <- as.vector(NULL)
        inter_col <- as.vector(unlist(target_data[["target"]][["knownDrugs"]][["rows"]][[y]][["drug"]]))
        inter_col <- as.character(inter_col)
        #make sure to only store unique records to 
        if(length(unique_chembl) == 0){
          unique_chembl <- inter_col[1]
          # #if drug has no approved indications, then add an N/A for the record
          # if(length(inter_col == 5)){
          #   inter_col[6] <- as.character("N/A")
          # }
          temp_row <- c(target_data[["target"]][["id"]],target_data[["target"]][["approvedSymbol"]],target_data[["target"]][["approvedName"]],inter_col,sig_dbs[i], sig_paths[i])
          temp_row <- as.data.frame(t(temp_row))
          merged_drugs <- as.data.frame(rbind(merged_drugs,temp_row),stringsAsFactors = FALSE, drop = FALSE)
          #write.table(df_init, file = outfile, row.names = FALSE, col.names=TRUE, sep = "\t", append = FALSE)
          
        }else if(grep(inter_col[1],unique_chembl) > 0){
          print(paste0("Duplicate drug found for same target: ",inter_col[1]))
          next()
        }else if(grep(inter_col[1],unique_chembl) == 0){
          unique_chembl <- c(unique_chembl, inter_col[1])
          #if drug has no approved indications, then add an N/A for the record
          if(length(inter_col == 5)){
            inter_col[6] <- as.character("N/A")
          }
          temp_row <- c(target_data[["target"]][["id"]],target_data[["target"]][["approvedSymbol"]],target_data[["target"]][["approvedName"]],inter_col,sig_dbs[i], sig_paths[i])
          temp_row <- as.data.frame(t(temp_row))
          merged_drugs <- as.data.frame(rbind(merged_drugs,temp_row),stringsAsFactors = FALSE)
        }
        
        rm(inter_col)
        rm(temp_row)
      }
      
      #$temp_df <- as.data.frame(do.call("cbind",list(target_data[["target"]][["id"]], target_data[["target"]][["approvedSymbol"]], target_data[["target"]][["approvedName"]], unlist(target_data[["target"]][["knownDrugs"]][["rows"]][[1]]),sig_dbs[i], sig_paths[i])))
      
      #drug_merged <- as.data.frame(do.call("cbind", list(target_data[["data"]][["target"]][["gene_info"]][["symbol"]],target_data[["data"]][["target"]][["gene_info"]][["name"]],target_data[["data"]][["target"]][["gene_info"]][["geneid"]],target_data[["data"]][["target"]][["activity"]],target_data[["data"]][["disease"]][["name"]],target_data[["data"]][["drug"]][["molecule_name"]],target_data[["data"]][["drug"]][["molecule_type"]],sig_dbs[i],sig_paths[i])))
      ##colnames(drug_merged) <- c("Symbol","Name","Ensembl_ID","Modulation","Disorder","Molecule_Name","Molecule_Type")
      #write.table(drug_merged, file = outfile, row.names = FALSE, col.names=FALSE, sep = "\t", append = TRUE)
    }
    rm(target_data)
  }
}
colnames(merged_drugs) <- c("Target_ID","Target_Symbol","Target_Name","Drug_ID","Drug_Name","Is_FDA_Approved","Highest_Clinical_Trial_Phase","Has_Been_Withdrawn","Pathway_DB","Pathway_Name")
merged_drugs <- unique(merged_drugs)
write.table(merged_drugs, file = outfile, row.names = FALSE, col.names=TRUE, sep = "\t", append = FALSE)
print("Done")  
