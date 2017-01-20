################################################################################################################
################################################################################################################
####
#### OMIM_Generation_Input_LOOCV: # Script that generates an input file for LOOCV process, from the data  
####                               gene-disease associations coming from OMIM. Downloaded via Biomart.
#### 
################################################################################################################
################################################################################################################
####
#### 1. NAME: OMIM_Generation_Input_LOOCV.R 
#### 2. CONTENTS: Protocol to generate an input file for the LOOCV process to be performed in other scripts
####              in order to compare the performance of our method with previous ones. 
#### 3. CREATED: 02/12/2016 by Alberto Valdeolivas.
####
################################################################################################################
################################################################################################################
#### 4. DESCRIPTION: Protocol to generate an input file for the LOOCV process to be performed in other scripts
####              in order to compare the performance of our method with previous ones. For the input we take
####              data of gene-disease association from OMIM. (Biomart). To achieve the LOOCV as defined in the 
####              manuscript we need to save just those diseases with at least two associated genes which are 
####              also present in our networks (So we can use these associated genes as seed nodes.)
####    4.1. INPUT PARAMETERS:
####        - Networks to be considered: We need to know the networks we would like to use to extract the 
####          proteins. We need seed genes to be present in our networks. Options: PPI: Protein protein
####          interaction network; PATH: Pathway network; COEX: Co-expression network. To be introduced 
####          separated by commas.
####        - Pool of nodes or Common nodes: Spicify if we want the ensemble of common nodes (Intersection) 
####          of different layers or the pool of them (Union).
####        - Output file Name: Two files will be generated. One data frame with Gene-Disease associations and
####          an another one containing the nodes that we are considering.
####    4.2. OMIM DATA: Downloaded from BIOMART (02/12/2016) Using biomaRt (biomaRt_2.30.0 Package verison)
################################################################################################################
################################################################################################################
#### 5. OUTPUT: A data frame containing two columns. In the first one, we have the Disgenet disease ID and in 
####            second one we have the genes associated with those diseases.
####    
################################################################################################################
################################################################################################################
################################################################################################################
#### 6. EXAMPLE OF EXECUTION: Rscript OMIM_Generation_Input_LOOCV.R PPI,PATH,COEX Union OMIM_All_Nodes
####                          
################################################################################################################
################################################################################################################

rm(list=ls())

################################################################################################################
## 1.- We load the common network generation functions, and we install/load the packages needed. 
################################################################################################################
source("LOOCV_Utils.R")
CRAN.packages <- c("igraph")
bioconductor.packages <- c("biomaRt")
install.packages.if.necessary(CRAN.packages, bioconductor.packages)

# We create a directoy called networks_files to save the files generated.
data.dir <- "OMIM_DATA"
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
## 2.3.- type of argument.
type <- args[2]
## 2.4.- Name for the output file.
Output_FileName_GeneDisease_txt <- paste("OMIM_DATA/InputFileLOOCV_",args[3],".txt",sep="",collapse = "")
Output_FileName_Nodes_txt <- paste("OMIM_DATA/",args[3],".txt",sep="",collapse = "")

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
## 4.- We get the data from BiomaRt and we just keep those interesting diseases.
################################################################################################################
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

Gene_Phenotype_relation <- getBM(attributes=c('mim_morbid','hgnc_symbol'), 
                                 filters='hgnc_symbol', values=Set_of_nodes, mart=ensembl)

### More than two. I take just the name of diseases with at least two.

OMIM_Atleast_2Genes_Set <- names(which(table(Gene_Phenotype_relation$mim_morbid) > 1))
OMIM_Atleast_2Genes_Set_df <- Gene_Phenotype_relation[which(Gene_Phenotype_relation$mim_morbid %in% OMIM_Atleast_2Genes_Set),]

## I order the data frame by disease.
OMIM_Atleast_2Genes_Set_df <- OMIM_Atleast_2Genes_Set_df[with(OMIM_Atleast_2Genes_Set_df, order(mim_morbid)),]

################################################################################################################
## 5.- We write the output files containing gene-disease associations and the nodes used as common.
################################################################################################################
write.table(OMIM_Atleast_2Genes_Set_df,Output_FileName_GeneDisease_txt,sep="\t",
            row.names = FALSE, dec=".", quote=FALSE)  

write.table(Set_of_nodes,Output_FileName_Nodes_txt,sep="\t",
            row.names = FALSE, dec=".", quote=FALSE,col.names = FALSE)  