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
outfile1 <- paste0(infile,"-RankedTargets.tsv")
setwd("~/fsl_groups/fslg_PickettLabGroup/spia")
#setwd("~/")
merged_drugs <- as.data.frame(NULL)
trac_vector <- as.vector(c(1,2,3,9,10,11,18,19,20,26,27,28))

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
  #i <- 7
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
    #print("in Reactome loop")
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
#  entrez.genes <- db_entrezID
#  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
#  genes <- getBM(
#    filters="entrezgene_id",
#    attributes=c("ensembl_gene_id", "entrezgene_id"),#, "hgnc_symbol", "uniprot_gn_id"),#, "go_id", "name_1006"),
#    values=db_entrezID,
#    mart=mart,
#    unique = TRUE)
  
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
      associatedDiseases {
        count
      }
      tractability {
        modality
        id
        value
      }
      safetyLiabilities {
        event
      }
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
        tractability <- 0
        unique_chembl <- as.vector(NULL)
        inter_col <- as.vector(unlist(target_data[["target"]][["knownDrugs"]][["rows"]][[y]][["drug"]]))
        inter_col <- as.character(inter_col)
        #sum tractability at "approved", "advanced clinical trials", or "phase 1 trials" for SM, Ab, PR, and OC
        for(z in 1:length(trac_vector)){
          #z <- 1
          if(length(target_data[["target"]][["tractability"]])>0){
              if(target_data[["target"]][["tractability"]][[z]][["value"]]=="TRUE"){
                tractability <- tractability+1
              }
          }
        }
        #make sure to only store unique records to 
        if(length(unique_chembl) == 0){
          unique_chembl <- inter_col[1]
          # #if drug has no approved indications, then add an N/A for the record
          # if(length(inter_col == 5)){
          #   inter_col[6] <- as.character("N/A")
          # }
          temp_row <- c(target_data[["target"]][["id"]],target_data[["target"]][["approvedSymbol"]],target_data[["target"]][["approvedName"]],target_data[["target"]][["associatedDiseases"]][["count"]],tractability,length(target_data[["target"]][["safetyLiabilities"]]),target_data[["target"]][["knownDrugs"]][["uniqueDrugs"]],inter_col,sig_dbs[i], sig_paths[i])
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
          temp_row <- c(target_data[["target"]][["id"]],target_data[["target"]][["approvedSymbol"]],target_data[["target"]][["approvedName"]],target_data[["target"]][["associatedDiseases"]][["count"]],tractability,length(target_data[["target"]][["safetyLiabilities"]]),target_data[["target"]][["knownDrugs"]][["uniqueDrugs"]],inter_col,sig_dbs[i], sig_paths[i])
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
#colnames(merged_drugs) <- c("Target_ID","Target_Symbol","Target_Name","Associated_Disease_Count","Tractability_Count","Safety_Liabilities","Number_Unique_Drugs","Drug_ID","Drug_Name","Is_FDA_Approved","Highest_Clinical_Trial_Phase","Has_Been_Withdrawn","Pathway_DB","Pathway_Name")
##colnames(merged_drugs) <- c("Target_ID","Target_Symbol","Target_Name","Drug_ID","Drug_Name","Is_FDA_Approved","Highest_Clinical_Trial_Phase","Has_Been_Withdrawn","Pathway_DB","Pathway_Name")
#merged_drugs <- unique(merged_drugs)
#write.table(merged_drugs, file = outfile, row.names = FALSE, col.names=TRUE, sep = "\t", append = FALSE)
print("Data gathering...complete")
print("Computing data...")


###########
#sort the existing table by the following order of criteria
# 1)  Pathway_Count,
# 2)  Tractability_Count,
# 2)  Number_Approved,
# 3)  Safety_Liabilities,
# 4)  Number_Unique_Drugs,
# 5)  Associated_Disease_Count,
# 6)  Number_phase3,
# 7)  Number_phase2,
# 8)  Number_phase1,
# 9)  Number_phase4

#Sorting step 1: iterate over whole table
#find counts for each target across all pathways
merged_drugs2 <- subset(merged_drugs, select = -c(Drug_ID,
                                                   Drug_Name,
                                                   Is_FDA_Approved,
                                                   Highest_Clinical_Trial_Phase,
                                                   Has_Been_Withdrawn
                                                   ))
merged_drugs2 <- unique(merged_drugs2)
merged_drugs2$Pathway_Name <- as.character(merged_drugs2$Pathway_Name)
merged_drugs2$Target_Symbol <- as.character(merged_drugs2$Target_Symbol)
target1 <- as.vector(merged_drugs2$Target_Symbol)
target1 <- unique(target1)
num_pathways <- as.vector(as.numeric(NULL))
for(a in 1:length(target1)){
  #a <- 3
  merged_drugs_temp <- as.data.frame(NULL)
  merged_drugs_temp <- subset(merged_drugs2, select = c(Pathway_Name,Target_Symbol))
  #merged_drugs_temp$Pathway_Name <- as.character(merged_drugs_temp$Pathway_Name)
  #merged_drugs_temp$Target_Symbol <- as.character(merged_drugs_temp$Target_Symbol)
  merged_drugs_temp <- subset(merged_drugs_temp, Target_Symbol == target1[a])
  num_pathways[a] <- nrow(merged_drugs_temp)
  #merged_drugs_temp <- subset(merged_drugs2, select = Pathway_Name & Target_Symbol %in% c(target1[a]))
}
rm(merged_drugs_temp)
rm(merged_drugs2)
num_path_df <- as.data.frame(cbind(num_pathways,target1))
num_path_df$num_pathways <- as.numeric(as.character(num_path_df$num_pathways))
num_path_df$target1 <- as.character(num_path_df$target1)
merged_drugs <- merge(merged_drugs,num_path_df, by.x='Target_Symbol', by.y='target1')

#Step2: iterate through each target
merged_drugs3 <- subset(merged_drugs, select = c(Target_Symbol,
                                                  Drug_Name,
                                                  Is_FDA_Approved,
                                                  Highest_Clinical_Trial_Phase,
                                                  Has_Been_Withdrawn
    ),
    stringsAsFactors = FALSE
  )
num4v <- as.vector(as.numeric(NULL))
num3v <- as.vector(as.numeric(NULL))
num2v <- as.vector(as.numeric(NULL))
num1v <- as.vector(as.numeric(NULL))
numApprovedv <- as.vector(as.numeric(NULL))

merged_drugs3 <- unique(merged_drugs3, stringsAsFactors = FALSE)
merged_drugs3$Target_Symbol<- as.character(merged_drugs3$Target_Symbol)
merged_drugs3$Target_Symbol<- as.character(merged_drugs3$Target_Symbol)
merged_drugs3$Is_FDA_Approved<- as.character(merged_drugs3$Is_FDA_Approved)
merged_drugs3$Highest_Clinical_Trial_Phase<- as.numeric(as.character(merged_drugs3$Highest_Clinical_Trial_Phase))

for(b in 1:length(target1)){
  #b <- 1
  phase_temp <- as.data.frame(subset(merged_drugs3, Target_Symbol == target1[b]), stringsAsFactors = FALSE)
  for(c in 1:nrow(phase_temp)){
    numApprovedv[b] <- as.numeric(sum(phase_temp$Is_FDA_Approved == "TRUE"))
    num3v[b] <- as.numeric(sum(phase_temp$Highest_Clinical_Trial_Phase == 3))
    num2v[b] <- as.numeric(sum(phase_temp$Highest_Clinical_Trial_Phase == 2))
    num1v[b] <- as.numeric(sum(phase_temp$Highest_Clinical_Trial_Phase == 1))
    num4v[b] <- as.numeric(sum(phase_temp$Highest_Clinical_Trial_Phase == 4))
  }
}
rm(merged_drugs3)
rm(phase_temp)
num_df <- as.data.frame(cbind(target1,numApprovedv, num3v, num2v, num1v, num4v))
merged_drugs <- merge(merged_drugs,num_df, by.x='Target_Symbol', by.y='target1')

colnames(merged_drugs) <- c("Target_Symbol","Target_ID","Target_Name","Associated_Disease_Count","Tractability_Count","Safety_Liabilities","Number_Unique_Drugs","Drug_ID","Drug_Name","Is_FDA_Approved","Highest_Clinical_Trial_Phase","Has_Been_Withdrawn","Pathway_DB","Pathway_Name","Target_in_Pathways","num_Approved_Drugs","num_Phase3","num_Phase2","num_Phase1","num_Phase4")
merged_drugs$Tractability_Count <- as.numeric(as.character(merged_drugs$Tractability_Count))
merged_drugs$Safety_Liabilities <- as.numeric(as.character(merged_drugs$Safety_Liabilities))
merged_drugs$Number_Unique_Drugs <- as.numeric(as.character(merged_drugs$Number_Unique_Drugs))
merged_drugs$num_Approved_Drugs <- as.numeric(as.character(merged_drugs$num_Approved_Drugs))
merged_drugs$num_Phase3 <- as.numeric(as.character(merged_drugs$num_Phase3))
merged_drugs$num_Phase2 <- as.numeric(as.character(merged_drugs$num_Phase2))
merged_drugs$num_Phase1 <- as.numeric(as.character(merged_drugs$num_Phase1))
merged_drugs$num_Phase4 <- as.numeric(as.character(merged_drugs$num_Phase4))
merged_drugs$Associated_Disease_Count <- as.numeric(as.character(merged_drugs$Associated_Disease_Count))

merged_drugs <- merged_drugs[order(-merged_drugs$Target_in_Pathways,
                                     -merged_drugs$Tractability_Count,
                                     -merged_drugs$num_Approved_Drugs,
                                     merged_drugs$Safety_Liabilities,
                                     -merged_drugs$Number_Unique_Drugs,
                                     -merged_drugs$Associated_Disease_Count,
                                     -merged_drugs$num_Phase3,
                                     -merged_drugs$num_Phase2,
                                     -merged_drugs$num_Phase1,
                                     -merged_drugs$num_Phase4
),]
write.table(merged_drugs, file = outfile, row.names = FALSE, col.names=TRUE, sep = "\t", append = FALSE)

merged_drugs1 <- as.data.frame(merged_drugs, stringsAsFactors = FALSE)
merged_drugs1 <- as.data.frame(select(merged_drugs1, select = -c("Drug_ID",
                                                                "Drug_Name",
                                                                "Is_FDA_Approved",
                                                                "Highest_Clinical_Trial_Phase",
                                                                "Has_Been_Withdrawn",
                                                                "Pathway_DB",
                                                                "Pathway_Name"
                                                                )
                                      ),stringsAsFactors = FALSE
                               )
merged_drugs1$Tractability_Count <- as.numeric(as.character(merged_drugs1$Tractability_Count))
merged_drugs1$Safety_Liabilities <- as.numeric(as.character(merged_drugs1$Safety_Liabilities))
merged_drugs1$Number_Unique_Drugs <- as.numeric(as.character(merged_drugs1$Number_Unique_Drugs))
merged_drugs1$num_Approved_Drugs <- as.numeric(as.character(merged_drugs1$num_Approved_Drugs))
merged_drugs1$num_Phase3 <- as.numeric(as.character(merged_drugs1$num_Phase3))
merged_drugs1$num_Phase2 <- as.numeric(as.character(merged_drugs1$num_Phase2))
merged_drugs1$num_Phase1 <- as.numeric(as.character(merged_drugs1$num_Phase1))
merged_drugs1$num_Phase4 <- as.numeric(as.character(merged_drugs1$num_Phase4))
merged_drugs1$Associated_Disease_Count <- as.numeric(as.character(merged_drugs1$Associated_Disease_Count))

merged_drugs1 <- unique(merged_drugs1)
merged_drugs1 <- merged_drugs1[order(-merged_drugs1$Target_in_Pathways,
                                     -merged_drugs1$Tractability_Count,
                                     -merged_drugs1$num_Approved_Drugs,
                                     merged_drugs1$Safety_Liabilities,
                                     -merged_drugs1$Number_Unique_Drugs,
                                     -merged_drugs1$Associated_Disease_Count,
                                     -merged_drugs1$num_Phase3,
                                     -merged_drugs1$num_Phase2,
                                     -merged_drugs1$num_Phase1,
                                     -merged_drugs1$num_Phase4
                                     ),]
write.table(merged_drugs1, file = outfile1, row.names = FALSE, col.names=TRUE, sep = "\t", append = FALSE)

print("Computing data...complete")
