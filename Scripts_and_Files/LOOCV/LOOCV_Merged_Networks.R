################################################################################################################
################################################################################################################
####
#### LOOCV_Merged_Metworks:  # Script that performs a validation of the ranking ability of RWR on Merged Networks 
####                           via a LOOCV process of disease-gene associations.
####                                          
################################################################################################################
################################################################################################################
####
#### 1. NAME: LOOCV_Merged_Metworks.R: 
#### 2. CONTENTS: Protocol to validate the ranking ability of RWR via a LOOCV when this method is performed
####              over monoplex networks merged on a multigraph (Graph with nodes that can be connected by
####              different type of nodes.)
#### 3. CREATED: 02/12/2016 by Alberto Valdeolivas.
####
################################################################################################################
################################################################################################################
#### 4. DESCRIPTION: Protocol to validate the ranking ability of RWR via a LOOCV when this method is performed
####              over a multigraph. It takes disease-gene association and for 
####              each disease leaves one gene out in turn. The remaining genes associated with the disease
####              are used as seeds to feed the RWR. Then, we save how the method ranked the removed gene. 
####              Every gene is ranked and that information is recorder in the output file. It will be used
####              in other scripts to plot the results via a cumulative distributive function (CDF) 
####
####    4.1. INPUT PARAMETERS:
####      - 1.- Networks to be considered: RWR can be performed over monoplex networks (PPI, PATH or COEX) or over
####        any multiplex network resulting from the combinitaion of those ones.
####      - 2.- The name of the file containing the nodes used as benchmark for performance comparisons. 
####        Reference nodes (For instace: Commo nodes resulting from intersection, pool of nodes resulting from union)
####      - 3.- The input file name. The file should contain diseases in the first column and genes its
####        associated genes in the second one. 
####      - 4.- The name of the output file that will contain each disease, its related genes and their rank.
####       
################################################################################################################
################################################################################################################
#### 5. OUTPUT: A data frame containing each disease tested along with its related genes and the rank for all
####            of them in RWR algorithm after performing the LOOCV. That file will be used to display the 
####            results.
################################################################################################################
################################################################################################################
#### 6. EXAMPLE OF EXECUTION: Rscript LOOCV_Merged_Networks.R PPI,PATH,COEX DISGENET_DATA/Pool_Nodes.txt DISGENET_DATA/InputFileLOOCV_Pool_Nodes.txt Results_32/Merged_PPI_PATH_COEX
####                          
################################################################################################################
################################################################################################################
rm(list=ls())

################################################################################################################
## 1.- We load the common network generation functions, and we install/load the packages needed. 
################################################################################################################
source("LOOCV_Utils.R")
CRAN.packages <- c("igraph","Matrix")
bioconductor.packages <- c()
install.packages.if.necessary(CRAN.packages, bioconductor.packages)
################################################################################################################
## 2.- We read the input arguments needed to run the script. 
################################################################################################################
print("Reading arguments...")
args <- commandArgs(trailingOnly = TRUE)

## 2.1.- Networks to be used.
# Network_List <- c("PPI","PATH","COEX")
Network_List <- unlist((strsplit(args[1],",")))

## 2.2.- Name of the file containing the reference nodes (Set of nodes to be considered)
# Set_of_nodes_fileName <- "DISGENET_DATA/Pool_Nodes.txt"
Set_of_nodes_fileName <- args[2]
Set_of_nodes <- read.csv(Set_of_nodes_fileName,stringsAsFactors = FALSE, header=FALSE)
Set_of_nodes <- Set_of_nodes$V1

## 2.3.- The name of the file containing the disease and their associated diseases.
# Input_FileName <- "DISGENET_DATA/InputFileLOOCV_Pool_Nodes.txt"
Input_FileName <- args[3]
Disease_Genes_Input <- read.csv(Input_FileName,sep="\t",stringsAsFactors = FALSE)
colnames(Disease_Genes_Input) <- c("Disease_ID","HGNC_Symbol")

## 2.4.- The Output file name. 
Output_FileName_txt <- paste(args[4],".txt",sep="",collapse = "")

rm(args)

################################################################################################################
## 3.- We merge the layers into a multigraph and we normalize the sparse matrix. 
################################################################################################################
print("Merging Networks...")
Merged_Network <- Merge.Layers(Network_List)
Merged_Adjacency_Matrix <- as_adjacency_matrix(Merged_Network,sparse=TRUE)

Total_Strengh_of_Nodes <- Matrix::colSums (Merged_Adjacency_Matrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)

print("Normalization...")
Merged_Adjacency_Matrix_Norm <- Normalize_Matrix_Sparse(Merged_Adjacency_Matrix,Total_Strengh_of_Nodes)


################################################################################################################
## 4.- We read all the records in the file containing diseases and their associated genes. For each disease
##     we leave one of its related genes out in turn. The rest are used as seed genes in the RWR. We look
##     for the ranking of the leave-out gene and we save that information. A file with the ranking of all the
##     genes is generated.
################################################################################################################# 

### We define a Data frame that will contain all the information.
LOOCV_Ranking_Results <- data.frame(Disease_ID=character(length = nrow(Disease_Genes_Input)),
                                    HGNC_Symbol = character(length = nrow(Disease_Genes_Input)),
                                    Global_Ranking=numeric(length = nrow(Disease_Genes_Input)),
                                    Total_Rankings = numeric(length = nrow(Disease_Genes_Input)),
                                    stringsAsFactors = FALSE)

### Parameter r. Global Restart probability of the RWR. It means that after every step we can continue with
### the walk or we can go back to the initial point. This value is always fixed to 0.7
############################
r <- 0.7

for(i in 1:nrow(Disease_Genes_Input)){
  
  ## We control where we are.
  print(paste("RWRs Percentage Completed: ", round(((i*100)/nrow(Disease_Genes_Input)),2),"%",sep="",collapse = NULL))
  
  ## We get the current Disease and Gene, and we take the genes that are associated to the same disease, 
  ## removing the current one. Those Genes will be used as 
  ## Seeds for our RWR on Multiplex. All of them will have the same restart probability parameter.
  Current_Disease <- Disease_Genes_Input$Disease_ID[i]
  Current_Gene <- Disease_Genes_Input$HGNC_Symbol[i]
  All_Disease_Genes <- Disease_Genes_Input$HGNC_Symbol[which(Disease_Genes_Input$Disease_ID %in% Current_Disease)]
  Seed_Genes <- All_Disease_Genes[which(!All_Disease_Genes %in% Current_Gene)]
  
  #### We prepare the Data frame with the seed genes and their scores to introduce to the RWR function
  SeedGenes <- data.frame(GeneNames = Seed_Genes, Score = 1/length(Seed_Genes),stringsAsFactors = FALSE)
  
  #### RWR - We apply the monoplex case with the Supra Matrix.
  Random_Walk_Results <- Random_Walk_Restart(Merged_Adjacency_Matrix_Norm,r,SeedGenes)
  
  #### We sort the results
  Global_results <- data.frame(GeneNames = row.names(Random_Walk_Results), Score = Random_Walk_Results[,1], stringsAsFactors = FALSE)
  Global_results <- Global_results[with(Global_results, order(-Score,GeneNames)), ]
  
  ### We remove the seed genes from the Ranking
  Global_results_NoSeeds <- Global_results[which(!Global_results$GeneNames %in% Seed_Genes),]
  
  ## Additionally we have to remove from the ranking the nodes that are not coming from the reference
  Global_results_NoSeeds_JustRefenceNodes <- Global_results_NoSeeds[which(Global_results_NoSeeds$GeneNames %in% Set_of_nodes),]
  
  ### We fill the data Frame with the results for each RWR
  LOOCV_Ranking_Results$Disease_ID[i] <- Current_Disease
  LOOCV_Ranking_Results$HGNC_Symbol[i] <- Current_Gene
  LOOCV_Ranking_Results$Global_Ranking[i] <- which(Global_results_NoSeeds_JustRefenceNodes$GeneNames ==Current_Gene)
  LOOCV_Ranking_Results$Total_Rankings[i] <- nrow(Global_results_NoSeeds_JustRefenceNodes)
}  

#### We write the results!
write.table(LOOCV_Ranking_Results,Output_FileName_txt,sep="\t",row.names = FALSE, dec=".",quote=FALSE)

  