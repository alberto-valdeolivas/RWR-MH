################################################################################################################
################################################################################################################
####
#### LOOCV_Multiplex_Heterogeneous:  # Script that performs a validation of the ranking ability of RWR 
####                                          via a LOOCV process of disease-gene associations. 
####                                          
################################################################################################################
################################################################################################################
####
#### 1. NAME: LOOCV_Multiplex_Heterogeneous.R: 
#### 2. CONTENTS: Protocol to validate the ranking ability of RWR via a LOOCV when this method is performed
####              over a multiplex and heterogeneous network. (It also works for monoplex heterogeneous networks.)
#### 3. CREATED: 05/12/2016 by Alberto Valdeolivas.
####
################################################################################################################
################################################################################################################
#### 4. DESCRIPTION: Protocol to validate the ranking ability of RWR via a LOOCV when this method is performed
####              over a multiplex and heterogeneous network. It takes disease-gene association and for 
####              each disease leaves one gene out in turn. The remaining genes associated with the disease
####              are used as seeds as well as the disease to feed the RWR. Then, we save how the method ranked 
####              the removed gene. Every gene is ranked and that information is recorder in the output file. 
####              It will be used in other scripts to plot the results via a cumulative distributive function (CDF) 
####              
####    4.1. INPUT PARAMETERS:
####      - 1.- Networks to be considered in the multiplex: RWR can be performed over monoplex networks (PPI, PATH or COEX) 
####        or over any multiplex network resulting from the combinitaion of those ones.
####      - 2.- The name of the file containing the nodes used as benchmark for performance comparisons. 
####        Reference nodes (For instace: Common nodes resulting from intersection, pool of nodes resulting from union)
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
#### 6. EXAMPLE OF EXECUTION: Rscript LOOCV_Multiplex_Heterogeneous.R PPI,PATH,COEX OMIM_DATA/PPI_Nodes.txt OMIM_DATA/InputFileLOOCV_PPI_Nodes.txt Results_32/Heterogeneous/PPI_Heterogeneous
####                          
################################################################################################################
################################################################################################################
rm(list=ls())

################################################################################################################
## 1.- We load the common network generation functions, and we install/load the packages needed. 
################################################################################################################
source("LOOCV_Utils.R")
CRAN.packages <- c("igraph","Rcpp","Matrix")
bioconductor.packages <- c("biomaRt")
install.packages.if.necessary(CRAN.packages, bioconductor.packages)
sourceCpp("Geometric_Mean.cpp")
# sourceCpp("get_transition_multiplex_rcpp.cpp")
################################################################################################################
## 2.- We read the input arguments needed to run the script. 
################################################################################################################
print("Reading arguments...")
args <- commandArgs(trailingOnly = TRUE)

## 2.1.- Networks to be used.
# Network_List <- c("PPI","PATH","COEX")
Network_List <- unlist((strsplit(args[1],",")))

## 2.2.- Name of the file containing the reference nodes (Set of nodes to be considered)
# Set_of_nodes_fileName <- "OMIM_DATA/OMIM_All_Nodes.txt"
Set_of_nodes_fileName <- args[2]
Set_of_nodes <- read.csv(Set_of_nodes_fileName,stringsAsFactors = FALSE, header=FALSE)
Set_of_nodes <- Set_of_nodes$V1

## 2.3.- The name of the file containing the disease and their associated diseases.
# Input_FileName <- "OMIM_DATA/InputFileLOOCV_OMIM_All_Nodes.txt"
Input_FileName <- args[3]
Disease_Genes_Input <- read.csv(Input_FileName,sep="\t",stringsAsFactors = FALSE)
colnames(Disease_Genes_Input) <- c("Disease_ID","HGNC_Symbol")

## 2.4.- The Output file name. 
Output_FileName_txt <- paste(args[4],".txt",sep="",collapse = "")

rm(args)

################################################################################################################
## 3.- We read the layers of our multiplex network. We generate its supra-adjacency matrix.
################################################################################################################
print("Reading Layers of the Multiplex Network...")

L <- length(Network_List)
List_Layers <- read.layers(Network_List)

### We get a pool of nodes. Nodes in any of the layers. Some nodes will belong to just some layers, while others will be on all of them.
pool_nodes <- pool.of.nodes(List_Layers)
pool_nodes_sorted <- sort(pool_nodes)

### We add to each layer the missing nodes with no connections (no edges for them)
List_Layers_Allnodes <- add.missing.nodes(List_Layers, pool_nodes)

### We have to check that all the layers have the same number of Vertex. We save N as the number of Nodes in all the layers(After adding those that were missing.)
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


################################################################################################################
## 4.- We generate the heterogeneous network.
################################################################################################################
print("Generating the Heterogeneous-Multiplex Network...")

print("Reading the disease-similarity network...")

## Reading Disease-Disease similarity Network
Disease_table <- read.table("../Network_Generation/networks_files/DiseaseSimilarity_2016-12-06.gr",sep=" ")
Disease_Network <- graph.data.frame(Disease_table,directed=FALSE)
Disease_Network <- igraph::simplify(Disease_Network, remove.multiple = TRUE, remove.loops = TRUE)
AdjMatrix_Diseases <- as_adjacency_matrix(Disease_Network,sparse = TRUE)
# M <- Size of the Disease-disease similarity graph (number of nodes)
M <- nrow(AdjMatrix_Diseases)
##  4.1.-We generate the bipartite graph
################################################################################################################
print("Generating the bipartite graph between genes and diseases...")
Gene_Phenotype_relation <- get.disease.gene.relations(pool_nodes)

Bipartite_matrix_and_report <- get.bipartite.graph(pool_nodes_sorted, colnames(AdjMatrix_Diseases), Gene_Phenotype_relation,N,M)
Bipartite_matrix <- Bipartite_matrix_and_report[[1]]
Error_log <- Bipartite_matrix_and_report[[2]]

## We expand the biparite graph to fit the multiplex dimensions.
## The biparti matrix has now (N x M)  dimensions. However we need it to have NL x M
## The genes in all the layers have to point to the diseases

print("Adapting the bipartite graph to the multiplex...")

SupraBipartiteMatrix <- expand.bipartite.graph(N,L,M,Bipartite_matrix)


##  4.2.-From the supra-adjacency matrix of the Multiplex system, the adjacency matrix of the disease-disease 
##       Similarity network and the expandend bipartite graph we calculate a transition matrix according to Li
##       and Patra (2010)
################################################################################################################
print("Transforming from adjacency matrices to transition matrix...")

### Parameter lambda. If a gene is connected with a disease in the bipartite graph, the walker can either hop 
### between subnetworks or stay in the same subgraph. To quantify this phenomena the parameter, called jumping probability, 
### was introduced. The larger the value of lambda the higher probability of changing between subnetworks. 
############################
lambda <- 0.5 

#### Transition Matrix for the inter-subnetworks links

Transition_Protein_Disease <- get.transition.protein.disease(N,L,M,SupraBipartiteMatrix,lambda)
Transition_Disease_Protein <- get.transition.disease.protein(N,L,M,SupraBipartiteMatrix,lambda)

#### Transition Matrix for the intra-subnetworks links
# Transition_Multiplex_Network <- get.transition.multiplex(N,L,lambda,SupraAdjacencyMatrix,SupraBipartiteMatrix)
Transition_Multiplex_Network <- get.transition.multiplex(N,L,lambda,Supra_Adj_Matrix,SupraBipartiteMatrix)
Transition_Disease_Network <- get.transition.disease(M,lambda,AdjMatrix_Diseases,SupraBipartiteMatrix)

### We generate the global transiction matrix and we return it.
Transition_Multiplex_Heterogeneous_Matrix_1 <- cbind(Transition_Multiplex_Network, Transition_Protein_Disease)
Transition_Multiplex_Heterogeneous_Matrix_2 <- cbind(Transition_Disease_Protein, Transition_Disease_Network)
Transition_Multiplex_Heterogeneous_Matrix <- rbind(Transition_Multiplex_Heterogeneous_Matrix_1,Transition_Multiplex_Heterogeneous_Matrix_2)

### Parameter eta.The parameter eta controls the probability of restarting in each subnetwork
############################
eta <- 0.5

################################################################################################################
## 5.- We read all the records in the file containing diseases and their associated genes. For each disease
##     we leave one of its related genes out in turn. The rest are used as seed genes in the RWR along with the  
##     the current disease. We look for the ranking of the leave-out gene and we save that information.
##     A file with the ranking of all the genes and diseases is generated.
################################################################################################################
 

### Parameter r. Global Restart probability of the RWR. It means that after every step we can continue with
### the walk or we can go back to the initial point. This value is always fixed to 0.7
############################
r <- 0.7

### Parameter tau. Restart probability within the layes. We can assume that the probability to restart in one layer
### is larger than the probability to restart in a different one. In this section we assume an equal restart probability
### Section Results 3.1. of the manuscript.
############################
tau <- rep(1/L,L)
# tau <- c(1/L,1.9/L,0.1/L)


## In this case we need to impose to perform the validation that the diseases (that we are also going to take as seed) are
## in our disease-disease similarity network. 
Disease_Genes_Input <- Disease_Genes_Input[which(Disease_Genes_Input$Disease_ID %in% colnames(AdjMatrix_Diseases)),]

### We define a Data frame that will contain all the information.
LOOCV_Ranking_Results <- data.frame(Disease_ID=character(length = nrow(Disease_Genes_Input)),
                                    HGNC_Symbol = character(length = nrow(Disease_Genes_Input)),
                                    Global_Ranking=numeric(length = nrow(Disease_Genes_Input)),
                                    Total_Rankings = numeric(length = nrow(Disease_Genes_Input)),
                                    stringsAsFactors = FALSE)


for(i in 1:nrow(Disease_Genes_Input)){
  
  ## We control where we are.
  print(paste("RWRs Percentage Completed: ", round(((i*100)/nrow(Disease_Genes_Input)),2),"%",sep="",collapse = NULL))
  
  ## We get the current Disease and Gene, and we take the genes that are associated to the same disease, removing the current one. Those Genes will be used as 
  ## Seeds for our RWR on Multiplex. All of them will have the same restart probability parameter.
  Current_Disease <- Disease_Genes_Input$Disease_ID[i]
  Current_Gene <- Disease_Genes_Input$HGNC_Symbol[i]
  All_Disease_Genes <- Disease_Genes_Input$HGNC_Symbol[which(Disease_Genes_Input$Disease_ID %in% Current_Disease)]
  Seed_Genes <- All_Disease_Genes[which(!All_Disease_Genes %in% Current_Gene)]
  Seed_Genes_Layer_Labeled <- character(length = length(Seed_Genes)*L)
  Seeds_Genes_Scores <- numeric(length = length(Seed_Genes)*L)
  
  ### I have to prepare the Seed Genes, including in the name of the genes the layer number.
  Current_Gene_Labeled <- character()
  for (k in 1:L){
    Current_Gene_Labeled <- c(Current_Gene_Labeled,paste(Current_Gene,k,sep="_",collapse = "") )
    for (j in 1:length(Seed_Genes)){
      Seed_Genes_Layer_Labeled[((k-1)*length(Seed_Genes))+ j] <- paste(Seed_Genes[j],k,sep="_",collapse = "") 
      Seeds_Genes_Scores[((k-1)*length(Seed_Genes))+ j] <- ((1-eta) * tau[k])/length(Seed_Genes)
    }
  }  
  
  ### In this case also the associated diseases are seeds in the sistem.
  
  # Seeds_Genes_Scores <- rep((1-eta) * (tau/length(Seed_Genes)), length(Seed_Genes))
  Seeds_Diseases_Scores <- eta
  All_Scores <- c(Seeds_Genes_Scores,Seeds_Diseases_Scores)
  
  #### We prepare the Data frame with the seed genes and their scores to introduce to the RWR function
  All_Seeds <- data.frame(SeedNames = c(Seed_Genes_Layer_Labeled,Current_Disease), Score = All_Scores,stringsAsFactors = FALSE)
  
  #### Additionally we have to remove the link in the biparti graph between the current gene and the curren disease.
  Heterogeneoux_Transitiction_Matrix_RemovedLink <- Transition_Multiplex_Heterogeneous_Matrix
  Heterogeneoux_Transitiction_Matrix_RemovedLink[Current_Gene_Labeled,as.character(Current_Disease)] <- 0
  Heterogeneoux_Transitiction_Matrix_RemovedLink[as.character(Current_Disease),Current_Gene_Labeled] <- 0
  
  
  #### RWR - We apply the monoplex case with the Supra Matrix.
  Random_Walk_Results <- Random_Walk_Restart(Heterogeneoux_Transitiction_Matrix_RemovedLink,r,All_Seeds)
  
  ## We remove from the results
  rank_global <- data.frame(GeneNames = character(length = N), Score = 0)
  rank_global$GeneNames <- gsub("_1", "", row.names(Random_Walk_Results)[1:N])
  
  ## We perform the geometric mean of the results in the differemt layers.
  rank_global$Score <- Geometric_Mean(as.vector(Random_Walk_Results[,1]),L,N)

  Global_results <- rank_global[with(rank_global, order(-Score, GeneNames)), ]
  
  ### We remove the seed genes from the Ranking
  Global_results_NoSeeds <- Global_results[which(!Global_results$GeneNames %in% Seed_Genes),]
  
  ## Additionally we have to remove from the ranking the nodes that are not coming from the referemce
  Global_results_NoSeeds_JustRefenceNodes <- Global_results_NoSeeds[which(Global_results_NoSeeds$GeneNames %in% Set_of_nodes),]
  
  ### We fill the data Frame with the results for each RWR
  LOOCV_Ranking_Results$Disease_ID[i] <- Current_Disease
  LOOCV_Ranking_Results$HGNC_Symbol[i] <- Current_Gene
  LOOCV_Ranking_Results$Global_Ranking[i] <- which(Global_results_NoSeeds_JustRefenceNodes$GeneNames ==Current_Gene)
  LOOCV_Ranking_Results$Total_Rankings[i] <- nrow(Global_results_NoSeeds_JustRefenceNodes)
}  

print(paste("Lambda: ", lambda))
print(paste("Delta: ", delta_parameter))
print(paste("eta: ", eta))
print(paste("restart Probability, r: ", r))
print(paste("tau: ", tau))

#### We write the results!
write.table(LOOCV_Ranking_Results,Output_FileName_txt,sep="\t",row.names = FALSE, dec=".",quote=FALSE)
