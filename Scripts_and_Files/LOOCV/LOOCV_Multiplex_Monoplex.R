################################################################################################################
################################################################################################################
####
#### LOOCV_Multiplex_Monoplex:  # Script that performs a validation of the ranking ability of RWR 
####                                          via a LOOCV process of disease-gene associations.
####                                          
################################################################################################################
################################################################################################################
####
#### 1. NAME: LOOCV_Multiplex_Monoplex.R: 
#### 2. CONTENTS: Protocol to validate the ranking ability of RWR via a LOOCV when this method is performed
####              over different monoplex or multiplex networks. 
#### 3. CREATED: 01/12/2016 by Alberto Valdeolivas.
####
################################################################################################################
################################################################################################################
#### 4. DESCRIPTION: Protocol to validate the ranking ability of RWR via a LOOCV when this method is performed
####              over different monoplex or multiplex networks. It takes disease-gene association and for 
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
#### 6. EXAMPLE OF EXECUTION: Rscript LOOCV_Multiplex_Monoplex.R PPI,PATH,COEX DISGENET_DATA/CommonNodes.txt DISGENET_DATA/InputFileLOOCV_CommonNodes.txt Results_31/Multiplex
####                          
################################################################################################################
################################################################################################################
rm(list=ls())

################################################################################################################
## 1.- We load the common network generation functions, and we install/load the packages needed. 
################################################################################################################
source("LOOCV_Utils.R")
CRAN.packages <- c("igraph","Rcpp","Matrix")
bioconductor.packages <- c()
install.packages.if.necessary(CRAN.packages, bioconductor.packages)
sourceCpp("Geometric_Mean.cpp")

################################################################################################################
## 2.- We read the input arguments needed to run the script. 
################################################################################################################
print("Reading arguments...")
args <- commandArgs(trailingOnly = TRUE)

## 2.1.- Networks to be used.
# Network_List <- c("PPI","PATH","COEX")
Network_List <- unlist((strsplit(args[1],",")))

## 2.2.- Name of the file containing the reference nodes (Set of nodes to be considered)
# Set_of_nodes_fileName <- "DISGENET_DATA/CommonNodes.txt"
Set_of_nodes_fileName <- args[2]
Set_of_nodes <- read.csv(Set_of_nodes_fileName,stringsAsFactors = FALSE, header=FALSE)
Set_of_nodes <- Set_of_nodes$V1

## 2.3.- The name of the file containing the disease and their associated diseases.
# Input_FileName <- "DISGENET_DATA/InputFileLOOCV_CommonNodes.txt"
Input_FileName <- args[3]
Disease_Genes_Input <- read.csv(Input_FileName,sep="\t",stringsAsFactors = FALSE)
colnames(Disease_Genes_Input) <- c("Disease_ID","HGNC_Symbol")

## 2.4.- The Output file name. 
Output_FileName_txt <- paste(args[4],".txt",sep="",collapse = "")

rm(args)

################################################################################################################
## 3.- We read the layers of our network and we generate the pool of nodes. Also the Supra-adjacency matrix.
##     We additionally define the elements and parameters needed for thr RWR on a multiplex network. 
################################################################################################################
L <- length(Network_List)
List_Layers <- read.layers(Network_List)

### We get a pool of nodes. Nodes in any of the layers.Some nodes belong to just one layer, others are part everywhere.
pool_nodes <- pool.of.nodes(List_Layers)

### We add to each layer the missing nodes with no connections (no edges for them). Isolated proteins. The definiton
### of a multiplex normally deals for layers containing the same nodes - that's why we proceed this way.
List_Layers_Allnodes <- add.missing.nodes(List_Layers, pool_nodes)

### We have to check that all the layers have the same number of Vertex. We save N as the number of Nodes in all the 
### layers(After adding those which were missing.)
N <- Get.Number.Nodes(List_Layers_Allnodes)


## 3.1.- We build the Transition matrix of the multiplex graph. it
##       accounts for the different possible transition of the partile in a multiplex graph.
################################################################################################################
### We generate the matrix that accounts for all the different possible transitions in a multiplex graphs.
### Basically: After and iteration of the algo, we can stay in the same layer (we will move to any of
### the neighbours of the current node) or we can move to the same current node but in a different layer.

### Parameter delta. Quantifies the probability to change among layers after an step of the algo. In this section
### we consider same probability to stay or to change. 
############################
delta_parameter <- 0.5


### Now we build our multiplex/monoplex transition Matrix and we normalize it by column (Supra_Adj_Matrix)

### We differentiate between the monoplex and the multiplex situation.
print("Generating the Multiplex Transition Matrix...")

if (L > 1){
  ### Multiplex.
  Supra_Adj_Matrix <- get.supra.adj.multiplex(List_Layers_Allnodes,delta_parameter,N)
} else{
  ### Monoplex.
  Supra_Adj_Matrix <- as_adjacency_matrix(List_Layers_Allnodes[[1]],sparse = TRUE)
  Supra_Adj_Matrix <- Supra_Adj_Matrix[order(rownames(Supra_Adj_Matrix)),order(colnames(Supra_Adj_Matrix))]
  colnames(Supra_Adj_Matrix) <- paste(colnames(Supra_Adj_Matrix),1,sep="_")
  rownames(Supra_Adj_Matrix) <- paste(rownames(Supra_Adj_Matrix),1,sep="_")
}
  
### We perform the column normalizaiton.
Total_Strengh_of_Nodes <-Matrix::colSums(Supra_Adj_Matrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)
print("Getting the Column Normalization...")
Supra_Adj_Matrix_Normalized <-  t(t(Supra_Adj_Matrix)/(Total_Strengh_of_Nodes))

################################################################################################################
## 4.- We read all the records in the file containing diseases and their associated genes. For each disease
##     we leave one of its related genes out in turn. The rest are used as seed genes in the RWR. We look
##     for the ranking of the leave-out gene and we save that information. A file with the ranking of all the
##     genes is generated.
################################################################################################################

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

### Parameter tau. Restart probability within the layes. We can assume that the probability to restart in one layer
### is larger than the probability to restart in a different one. In this section we assume an equal restart probability
### Section Results 3.1. of the manuscript.
############################
tau <- rep(1/L,L)

for(i in 1:nrow(Disease_Genes_Input)){
  
  ## We print the progress of the process.
  print(paste("RWRs Percentage Completed: ", round(((i*100)/nrow(Disease_Genes_Input)),2),"%",sep="",collapse = NULL))
  
  ## We get the current Disease and Gene, and we take the genes that are associated to the same disease, removing the current one. 
  ## Those Genes will be used as Seeds for our RWR on Multiplex. All of them will have the same restart probability parameter.
  Current_Disease <- Disease_Genes_Input$Disease_ID[i]
  Current_Gene <- Disease_Genes_Input$HGNC_Symbol[i]
  All_Disease_Genes <- Disease_Genes_Input$HGNC_Symbol[which(Disease_Genes_Input$Disease_ID %in% Current_Disease)]
  Seed_Genes <- All_Disease_Genes[which(!All_Disease_Genes %in% Current_Gene)]
  Seed_Genes_Layer_Labeled <- character(length = length(Seed_Genes)*L)
  Seed_Genes_Layer_Score <- numeric(length = length(Seed_Genes)*L)
  
  ### I have to prepare the Seed Genes, including in the name of the genes the layer number.
  for (k in 1:L){
    for (j in 1:length(Seed_Genes)){
      Seed_Genes_Layer_Labeled[((k-1)*length(Seed_Genes))+ j] <- paste(Seed_Genes[j],k,sep="_",collapse = "")
      Seed_Genes_Layer_Score[((k-1)*length(Seed_Genes))+ j] <- tau[k]/(length(Seed_Genes)*L)
    }
  }  
  
  #### We prepare the Data frame with the seed genes and their scores to introduce to the RWR function
  SeedGenes <- data.frame(GeneNames = Seed_Genes_Layer_Labeled, Score = Seed_Genes_Layer_Score,stringsAsFactors = FALSE)
  
  #### RWR - We apply the monoplex case with the Supra Adjacency Matrix.
  Random_Walk_Results <- Random_Walk_Restart(Supra_Adj_Matrix_Normalized,r,SeedGenes)
  
  ## To generate the ranking we just take the genes in the first layer; Then we calculate 
  ## the geometric mean taking in account results for each gene in every layer.
  rank_global <- data.frame(GeneNames = character(length = N), Score = 0)
  rank_global$GeneNames <- gsub("_1", "", row.names(Random_Walk_Results)[1:N])
  rank_global$Score <- Geometric_Mean(as.vector(Random_Walk_Results[,1]),L,N)
  
  Global_results <- rank_global[with(rank_global, order(-Score, GeneNames)), ]
  
  
  ### We remove the seed genes from the Ranking
  Global_results_NoSeeds <- Global_results[which(!Global_results$GeneNames %in% Seed_Genes),]
  
  ## Additionally we have to remove from the ranking the nodes that are not coming from the reference nodes.
  Global_results_NoSeeds_JustRefenceNodes <- Global_results_NoSeeds[which(Global_results_NoSeeds$GeneNames %in% Set_of_nodes),]
  
  
  ### We fill the data Frame with the results for each RWR
  LOOCV_Ranking_Results$Disease_ID[i] <- Current_Disease
  LOOCV_Ranking_Results$HGNC_Symbol[i] <- Current_Gene
  LOOCV_Ranking_Results$Global_Ranking[i] <- which(Global_results_NoSeeds_JustRefenceNodes$GeneNames ==Current_Gene)
  LOOCV_Ranking_Results$Total_Rankings[i] <- nrow(Global_results_NoSeeds_JustRefenceNodes)
}  

#### We write the results!
write.table(LOOCV_Ranking_Results,Output_FileName_txt,sep="\t",row.names = FALSE, dec=".",quote=FALSE)
