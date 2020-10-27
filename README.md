# Pathway2Targets

## ARMOR
Instructions on downloading and installing ARMOR can be found here: https://github.com/csoneson/ARMOR.

## Converting Ensembl IDs to Entrez IDs
After successfully completing a run with ARMOR, the edgeR output file should be used as input to the Get_entrezID_from_ENSG.R. This script retrieves the relevant NCBI Entrez Gene IDs for the Ensembl Gene IDs output by ARMOR. 

## Running SPIA
The output of the Get_entrezID_from_ENSG.R script can then be used as input for the Signaling Pathway Impact Analysis (SPIA) algorithm using the SPIA_Code.Rmd script.

## Running Pathway2Targets
The output of the SPIA_Code.Rmd script should be used as input to the Pathways2Targets.R script, which will identify nodes in each significant pathway and use the opentargets.org database to identify nodes that are known drug targets.
