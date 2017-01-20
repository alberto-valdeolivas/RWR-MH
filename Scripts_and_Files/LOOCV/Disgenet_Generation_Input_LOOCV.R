################################################################################################################
################################################################################################################
####
#### Disgenet_Generation_Input_LOOCV: # Script that generates an input file for LOOCV process, from the data  
####                               gene-disease associations coming from DISGENET.
#### 
################################################################################################################
################################################################################################################
####
#### 1. NAME: Disgenet_Generation_Input_LOOCV.R 
#### 2. CONTENTS: Protocol to generate an input file for the LOOCV process to be performed in other scripts
####              in order to compare the performance of our method with previous ones. 
#### 3. CREATED: 01/12/2016 by Alberto Valdeolivas.
####
################################################################################################################
################################################################################################################
#### 4. DESCRIPTION: Protocol to generate an input file for the LOOCV process to be performed in other scripts
####              in order to compare the performance of our method with previous ones. For the input we take
####              data of gene-disease association from DISGENET. To achieve the LOOCV as defined in the 
####              manuscript we need to save just those diseases with at least two associated genes which are 
####              also present in our networks (So we can use these associated genes as seed nodes.)
####    4.1. INPUT PARAMETERS:
####        - Networks to be considered: We need to know the networks we would like to use to extract the 
####          proteins. We need seed genes to be present in our networks. Options: PPI: Protein protein
####          interaction network; PATH: Pathway network; COEX: Co-expression network. To be introduced 
####          separated by commas.
####        - Score cutoff for DISGENET: Each gene-disease association in DISGENET has a condifent score. We 
####          need to introduce this parameter. (We recommend 0.15)
####        - Pool of nodes or Common nodes: Spicify if we want the ensemble of common nodes (Intersection) 
####          of different layers or the pool of them (Union).
####        - Output file Name: Two files will be generated. One data frame with Gene-Disease associations and
####          an another one containing the nodes that we are considering.
####    4.2. DISGENET DATA:
####        - DOWNLOADED FROM: (01/12/2016). 
####        "http://www.disgenet.org/ds/DisGeNET/results/curated_gene_disease_associations.tsv.gz"
####        This file corresponds to DisGeNET v4.0 (created on Jun 07 2016)
################################################################################################################
################################################################################################################
#### 5. OUTPUT: A data frame containing two columns. In the first one, we have the Disgenet disease ID and in 
####            second one we have the genes associated with those diseases.

################ COMMENTS: 
### Genes from disgenet. // ## Info about the file version and downloaded day in the top comments of the file.
### Score conversation. We used a Score cutoff of 0.15 in this version, but it is a parameter so we can modify it
### https://thinklab.com/discussion/processing-disgenet-for-disease-gene-relationships/105
################################################################################################################
################################################################################################################
#### 6. EXAMPLE OF EXECUTION: Rscript Disgenet_Generation_Input_LOOCV.R PPI,PATH,COEX 0.15 Intersection CommonNodes
####                          
################################################################################################################
################################################################################################################

rm(list=ls())

################################################################################################################
## 1.- We load the common network generation functions, and we install/load the packages needed. 
################################################################################################################
source("LOOCV_Utils.R")
CRAN.packages <- c("igraph")
bioconductor.packages <- c()
install.packages.if.necessary(CRAN.packages, bioconductor.packages)

# We create a directoy called networks_files to save the files generated.
data.dir <- "DISGENET_DATA"
if (!file.exists(data.dir)) {
  dir.create(data.dir)
}

################################################################################################################
## 2.- We read the input arguments.
################################################################################################################
print("Reading arguments...")
args <- commandArgs(trailingOnly = TRUE)

## 2.1.- Networks to be used.
# Network_List <- c("PPI","PATH","COEX")
Network_List <- unlist((strsplit(args[1],",")))
## 2.2.- Score Cutoff.
# Score_cutoff <- 0.15
Score_cutoff <- as.numeric(args[2])
## 2.3.- type of argument.
type <- args[3]
## 2.4.- Name for the output file.
Output_FileName_GeneDisease_txt <- paste("DISGENET_DATA/InputFileLOOCV_",args[4],".txt",sep="",collapse = "")
Output_FileName_Nodes_txt <- paste("DISGENET_DATA/",args[4],".txt",sep="",collapse = "")

rm(args)
################################################################################################################
## 3.- We read the layers of the multiplex to generate the union or intersection of nodes.
################################################################################################################

List_Layers <- read.layers(Network_List)

if (type == "Intersection"){
  Set_of_nodes <- common.nodes(List_Layers)
} else {
  if (type == "Union") {
    Set_of_nodes <- pool.of.nodes(List_Layers)  
  } else {
    stop("Not correct Type Name")
  }
}

################################################################################################################
## 4.- We read the disgenet input file. Info about the version on the top of the file.
################################################################################################################
DisgeNet_Gen_Disease <- read.csv(file = "DISGENET_DATA/curated_gene_disease_associations.tsv", comment.char = "#",
                                 sep = "\t", stringsAsFactors = FALSE)

DisgeNet_Gen_Disease_score <- DisgeNet_Gen_Disease[which(DisgeNet_Gen_Disease$score >= Score_cutoff),]

DisgeNet_Gen_Disease_Set <- DisgeNet_Gen_Disease_score[which(DisgeNet_Gen_Disease_score$geneName %in% Set_of_nodes),]

### More than two. I take just the name of diseases with at least two.

Diseases_Atleast_2Genes_Set <- names(which(table(DisgeNet_Gen_Disease_Set$diseaseId) > 1))
Diseases_Atleast_2Genes_Set_df <- DisgeNet_Gen_Disease_Set[which(DisgeNet_Gen_Disease_Set$diseaseId %in% Diseases_Atleast_2Genes_Set),]

## I just keep the Disease_Gene Associations.
DisgeNet_Disease_Gen_df <- data.frame(Disease_ID =Diseases_Atleast_2Genes_Set_df$diseaseId, Gene_ID = Diseases_Atleast_2Genes_Set_df$geneName,
                                      stringsAsFactors = FALSE)

################################################################################################################
## 5.- We write the output files containing gene-disease associations and the nodes used as common.
################################################################################################################
write.table(DisgeNet_Disease_Gen_df,Output_FileName_GeneDisease_txt,sep="\t",
            row.names = FALSE, dec=".", quote=FALSE)  

write.table(Set_of_nodes,Output_FileName_Nodes_txt,sep="\t",
            row.names = FALSE, dec=".", quote=FALSE,col.names = FALSE)  


