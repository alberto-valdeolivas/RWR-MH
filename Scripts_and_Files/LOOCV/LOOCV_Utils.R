################################################################################################################
################################################################################################################
####
#### LOOCV_Utils: COMMON FUNCTIONS/APPLICATIONS USED TO GENERATE TO PERFORM THE LEAVE-ONE-OUT CROSS VALIDATION.
####              ALSO THE FUNCTIONS USED TO GENERATE THE INPUT FILES FOR THESE VALIDATIONS ARE PRESENT HERE.
#### 
################################################################################################################
################################################################################################################

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# GENERAL APPLICATION FUNCTIONS
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 1.- Package installation and load of libraries.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

install.packages.if.necessary <- function(CRAN.packages=c(), bioconductor.packages=c()) {
  for (p in CRAN.packages) {	
    if (!require(p, character.only=T)) { 	
      install.packages(p) 	
      library(p, character.only=T)  	
    }	
  }
  if (length(bioconductor.packages) > 0) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install()
  }
  for (p in bioconductor.packages) {
    if (!require(p, character.only=T)) {
      BiocManager::install(p)
      library(p, character.only=T)
    }
  }
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 2.- RANDOM WALK WITH RESTART
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

Random_Walk_Restart <- function(Network_Matrix, r,SeedGenes ){
  
  ### We define the threshold and the number maximum of iterations for the randon walker.
  Threeshold <- 1e-10
  NetworkSize <- ncol(Network_Matrix)
  
  ### We initialize the variables to control the flux in the RW algo.
  residue <- 1
  iter <- 1
  
  #### We define the prox_vector(The vector we will move after the first RW iteration. We start from The seed. We have to take in account
  #### that the walker with restart in some of the Seed genes, depending on the score we gave in that file).
  prox_vector <- matrix(0,nrow = NetworkSize,ncol=1)
  
  prox_vector[which(colnames(Network_Matrix) %in% SeedGenes[,1])] <- (SeedGenes[,2])
  
  prox_vector  <- prox_vector/sum(prox_vector)
  restart_vector <-  prox_vector
  
  while(residue >= Threeshold){
    
    old_prox_vector <- prox_vector
    prox_vector <- (1-r)*(Network_Matrix %*% prox_vector) + r*restart_vector
    residue <- sqrt(sum((prox_vector-old_prox_vector)^2))
    iter <- iter + 1; 
  }
  return(prox_vector) 
} 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# MULTIPLEX RELATED FUNCTIONS
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 3.- We read the different layers that integrate the multiplex network.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

read.layers <- function(vec_layers){
  
  ## "pre-allocate" an empty list of length size (Number of Layers)
  size <- length(vec_layers)
  Layers <- vector("list", size)
  
  ## We read the different Networks (Layers) of our Multiplex network. We also simplify the networks
  ## by removing possible self loops and possible multiple nodes.
  
  for (i in 1:size){
    if (vec_layers[i]=="PPI"){
      
      PPI_table <- read.table("../Network_Generation/networks_files/PPI_2016-11-23.gr",sep=" ")
      PPI_Network <- graph.data.frame(PPI_table,directed=FALSE)
      PPI_Network <- igraph::simplify(PPI_Network, remove.multiple = TRUE, remove.loops = TRUE)
      
      Layers[[i]] <- PPI_Network
      names(Layers)[i] <- "PPI_NETWORK"
    } else {
      if (vec_layers[i]=="PATH"){
        
        Pathway_table <- read.table("../Network_Generation/networks_files/AllPathways_2016-11-24.gr",sep=" ")
        Pathway_Network <- graph.data.frame(Pathway_table,directed=FALSE)
        Pathway_Network <- igraph::simplify(Pathway_Network, remove.multiple = TRUE, remove.loops = TRUE)
        
        Layers[[i]] <- Pathway_Network
        names(Layers)[i] <- "PATHWAY_NETWORK"
        
     } else {
          if (vec_layers[i]=="COEX"){
            Coex_table <- read.table("../Network_Generation/networks_files/Co-Expression_2016-11-23.gr",sep=" ")
            Coex_Network <- graph.data.frame(Coex_table,directed=FALSE)
            Coex_Network <- igraph::simplify(Coex_Network, remove.multiple = TRUE, remove.loops = TRUE)
            
            Layers[[i]] <- Coex_Network
            names(Layers)[i] <- "COEXPRESION_NETWORK"
          } else {
            stop("Not correct Network Name")
          }
      } 
    }
  }
  ## We return a list containing an igraph object for each layer.
  return(Layers)
} 

### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 4.- We obtain the common nodes of the considered layers of the multiplex.
#     Intersection of nodes that are present in all the considered layers.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

common.nodes <- function(Layers){
  List_nodes <- list()
  
  for (i in 1:length(Layers)){
    List_nodes[[i]] <- V(Layers[[i]])$name
    names(List_nodes)[i] <- names(Layers)[i]
    if (i == 1){
      Common_Nodes<- List_nodes[[i]] 
    } else {
      Common_Nodes <- Common_Nodes[which(Common_Nodes %in% List_nodes[[i]])]
    }
  }
  return(Common_Nodes)
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 5.- We obtain the pool of nodes of the considered layers of the multiplex.
#     (Union of nodes of all the layers)
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

pool.of.nodes <- function(Layers){
  
  ## We get the number of layers
  Nr_Layers <- length(Layers)
  
  ## We get the nodes of all the layers of the multiplex network. We save them into a vector.
  Node_Names_all <- character()
  for (i in 1:Nr_Layers) {
    Node_Names_Layer <- V(Layers[[i]])$name
    Node_Names_all <-c(Node_Names_all,Node_Names_Layer)
  }
  
  ## We remove duplicates.
  Node_Names_all <- unique(Node_Names_all)
  
  return(Node_Names_all)
} 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 6.- From the pool of nodes we add the missing proteins to each layer as 
#     isolated nodes.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
add.missing.nodes <- function (Layers,NodeNames) {
  
  ## We get the number of layers
  Nr_Layers <- length(Layers)
  
  ## We generate a new list of layers.
  Layers_New <- vector("list", Nr_Layers)
  
  ## We add to each layer the missing nodes of the total set of nodes, of the pool of nodes.
  for (i in 1:Nr_Layers){
    Node_Names_Layer <- V(Layers[[i]])$name
    Missing_Nodes <- NodeNames[which(!NodeNames %in% Node_Names_Layer)]
    Layers_New[[i]] <- add_vertices(Layers[[i]] ,length(Missing_Nodes), name=Missing_Nodes)
  }
  return(Layers_New)
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 7.- We check the total number of nodes in every layer and we return it. 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
Get.Number.Nodes <- function(Layers_Allnodes) {
  
  ## We get the number of layers
  Nr_Layers <- length(Layers_Allnodes)
  vector_check <- numeric(length = Nr_Layers)  
  
  for (i in 1:Nr_Layers){
    vector_check[i] <- vcount(Layers_Allnodes[[i]])  
  }
  
  if (all(vector_check == vector_check[1])){
    print("Number of nodes in every layer updated...")
    return(vector_check[1])  
  } else {
    stop("Not correct number of nodes in each Layer...")
  }
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 8.- We check the total number of nodes in every layer and we return it. 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

get.supra.adj.multiplex <- function(Layers,delta,N){
  
  ## IDEM_MATRIX.
  Idem_Matrix <- Diagonal(N, x = 1)
  L <- length(Layers)
  
  SupraAdjacencyMatrix <- Matrix(0,ncol=N*L,nrow=N*L,sparse = TRUE)
  
  Col_Node_Names <- character()
  Row_Node_Names <- character()
  
  ## We differentiate here between the monoplex situation and the multiplex one since delta parameter does not make sense as defined for a monoplex.
  
  for (i in 1:L){
    Adjacency_Layer <-  as_adjacency_matrix(Layers[[i]],sparse = TRUE)
    
    ## We order the matrix by the node name. This way all the matrix will have the same. Additionally we include a label with the layer number for each node name.
    Adjacency_Layer <- Adjacency_Layer[order(rownames(Adjacency_Layer)),order(colnames(Adjacency_Layer))]
    Layer_Col_Names <- paste(colnames(Adjacency_Layer),i,sep="_")
    Layer_Row_Names <- paste(rownames(Adjacency_Layer),i,sep="_")
    Col_Node_Names <- c(Col_Node_Names,Layer_Col_Names)
    Row_Node_Names <- c(Row_Node_Names,Layer_Row_Names)
    
    
    ## We fill the diagonal blocks with the transition probability within a layer
    Position_ini_row <- 1 + (i-1)*N
    Position_end_row <- N + (i-1)*N
    SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_row:Position_end_row)] <- (1-delta)*(Adjacency_Layer)
    
    ## We fill the off-diagonal blocks with the transition probability among layers.
    for (j in 1:L){
      Position_ini_col <- 1 + (j-1)*N
      Position_end_col <- N + (j-1)*N
      if (j != i){
        SupraAdjacencyMatrix[(Position_ini_row:Position_end_row),(Position_ini_col:Position_end_col)] <- (delta/(L-1))*Idem_Matrix
      }
    }
  } 
  
  
  rownames(SupraAdjacencyMatrix) <- Row_Node_Names
  colnames(SupraAdjacencyMatrix) <- Col_Node_Names
  
  return(SupraAdjacencyMatrix)
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# MULTIGRAPH RELATED FUNCTIONS
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 9.- Function that takes a vector with names of the layers of our multiplex 
#     network and generate their combined(merged) 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####


Merge.Layers <- function(vec_layers){
  
  ## We read the different Networks (Layers) of our Multiplex network.
  size <- length(vec_layers)
  Merged_table <- data.frame()
  
  for (i in 1:size){
    if (vec_layers[i]=="PPI"){
      
      PPI_table <- read.table("../Network_Generation/networks_files/PPI_2016-11-23.gr",sep=" ")
      Merged_table <- rbind(Merged_table,PPI_table)    
    } else {
      if (vec_layers[i]=="PATH"){
        
        Pathway_table <- read.table("../Network_Generation/networks_files/AllPathways_2016-11-24.gr",sep=" ")
        Merged_table <- rbind(Merged_table,Pathway_table) 
        
      } else {
        if (vec_layers[i]=="COEX"){
          Coex_table <- read.table("../Network_Generation/networks_files/Co-Expression_2016-11-23.gr",sep=" ")
          Merged_table <- rbind(Merged_table,Coex_table) 
            
        } else {
          stop("Not correct Network Name")
        }
      } 
    }
  }
  
  Merged_Network <- graph.data.frame(Merged_table,directed=FALSE)
  Merged_Network <- igraph::simplify(Merged_Network, remove.multiple = FALSE, remove.loops = TRUE)
  
  
  return(Merged_Network)
} 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 10.- Function that normalizes the weigthed sparse matrix(Merged Networks) 
#     This version is for Sparse matrix to save space in RAM memory when executing.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
Normalize_Matrix_Sparse <- function(Network_Matrix,Strengh_of_Nodes){
  
  Normalized_Weighted_Adj_Matrix <- Matrix(0,nrow = nrow(Network_Matrix), ncol = ncol(Network_Matrix),
                                    dimnames=list(colnames(Network_Matrix),rownames(Network_Matrix)),sparse = TRUE)
  for (i in 1:ncol(Network_Matrix)){
    Normalized_Weighted_Adj_Matrix[,i] <- (Network_Matrix[,i])/(Strengh_of_Nodes[i])
  }
  return(Normalized_Weighted_Adj_Matrix)
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# HETEROGENEOUS NETWORK RELATED FUNCTIONS
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 11.- Function that recovers the gene-phenotype relations for a set of proteins
#     retriving their associates diseases from Biomart (mim codes.)
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
get.disease.gene.relations <- function(proteins){
  
  print("Getting from Ensamble Gene-disease relations...")
  # We do not retrieve the data from BiomaRt at each execution in order to fix, the input.
  # We downloaded Gene-Phenotype associations for our pool of nodes of our multiplex network.
  # With that we can match with the proteins in the current execution. 
  
  # ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
  
  # Gene_Phenotype_relation <- getBM(attributes=c('hgnc_symbol','mim_morbid'), 
  #                                 filters='hgnc_symbol', values=proteins, mart=ensembl)
  # write.table(Gene_Phenotype_relation,"../Network_Generation/networks_files/Gene_Phenotype_relation.txt", 
  #            quote=FALSE, sep="\t",row.names = FALSE, col.names = TRUE)
  
  Gene_Phenotype_relation <- read.table("../Network_Generation/networks_files/Gene_Phenotype_relation.txt", 
                                        sep="\t", header=TRUE,stringsAsFactors = FALSE)
  
  ## We remove Gene-phenotype relations that are not in our proteins file. 
  Gene_Phenotype_relation <- Gene_Phenotype_relation[which(Gene_Phenotype_relation$hgnc_symbol %in% proteins), ]
  
  return(Gene_Phenotype_relation)
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 12.- Function that generates a bipartite graph from gene-disease relations.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####

get.bipartite.graph <- function(proteins_sorted, disease_names, Gene_Phenotype_relation,Number_Proteins,Number_Diseases){
  Bipartite_matrix <- Matrix(data=0, nrow=Number_Proteins, ncol=Number_Diseases)
  rownames(Bipartite_matrix) <- proteins_sorted
  colnames(Bipartite_matrix) <- disease_names
  log <- character()
  
  for (i in 1:Number_Proteins){
    current_gene <- proteins_sorted[i]
    current_mim <- Gene_Phenotype_relation$mim_morbid[which(Gene_Phenotype_relation$hgnc_symbol == current_gene)]
    
    # We have to check if the gene is the Gene_Phenotype_relation downloaded from Biomart. (Some weird HGNC symbols, not retrieved.)
    if (length(current_mim) > 0){
      for (j in 1:length(current_mim)){
        if (!is.na(current_mim[j])){
          # We need to identify the phenotypes position on the matrix.
          index_disease <- which(colnames(Bipartite_matrix) %in%  current_mim[j])
          # We have to check if that index is present in the matrix.
          if (length(index_disease) == 1){ 
            Bipartite_matrix[i,index_disease] <- 1
          } else {
            error_message <- paste("MIM_CODE", current_mim[j], length(index_disease), "No phenotype found",sep=";", collapse = NULL)
            log <- c(log,error_message)
          }           
        }
      }  
    } else {
      error_message <- paste("HGNC_Symbol", current_gene, "No HGNC found in Biomart",sep=";", collapse = NULL)
      log <- c(log,error_message)
    }
  }
  Bipartite_and_errorlog <- list(Bipartite_matrix,log)
  return(Bipartite_and_errorlog)
}

################################################################################################################
# 13.- We expand the bipartite graph to fit the dimensions of our multilpex system. From every layer, we can
# jump to the other subnetwork (disease similarity network in this case)
################################################################################################################
expand.bipartite.graph <- function(Number_Proteins,Number_Layers,Number_Diseases,Bipartite_matrix){
  
  SupraBipartiteMatrix <- Matrix(0,nrow=Number_Proteins*Number_Layers,ncol=Number_Diseases,sparse = TRUE)
  Row_Node_Names <- character()
  
  for (i in 1:Number_Layers){
    Layer_Row_Names <- paste(rownames(Bipartite_matrix),i,sep="_")
    Row_Node_Names <- c(Row_Node_Names,Layer_Row_Names)
    Position_ini_row <- 1 + (i-1)*Number_Proteins
    Position_end_row <- Number_Proteins + (i-1)*Number_Proteins
    SupraBipartiteMatrix[(Position_ini_row:Position_end_row),] <- Bipartite_matrix
  }  
  
  rownames(SupraBipartiteMatrix) <- Row_Node_Names
  colnames(SupraBipartiteMatrix) <- colnames(Bipartite_matrix)
  return(SupraBipartiteMatrix)
}

################################################################################################################
# 14.- We compute the transition matrices for the multiplex-heterogeneous system. 
#      Those matrices will generate the final transition matrices. 
################################################################################################################

# 14.1.-Protein-Disease Transition Matrix.
################################################################################################################
get.transition.protein.disease <- function(Number_Proteins,Number_Layers,Number_Diseases,SupraBipartiteMatrix,lambda){
  
  Transition_Protein_Disease <- Matrix(0,nrow=Number_Proteins*Number_Layers,ncol=Number_Diseases,sparse = TRUE)
  colnames(Transition_Protein_Disease) <- colnames(SupraBipartiteMatrix)
  rownames(Transition_Protein_Disease) <- rownames(SupraBipartiteMatrix)
  
  Col_Sum_Bipartite <- colSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)
  
  for (j in 1:Number_Diseases){
    if (Col_Sum_Bipartite[j] != 0){
      Transition_Protein_Disease[,j] <- (lambda*SupraBipartiteMatrix[,j]) /Col_Sum_Bipartite[j]
    }
  }
  return(Transition_Protein_Disease)
}

# 14.2.-Disease-Protein Transition Matrix.
################################################################################################################
get.transition.disease.protein <- function(Number_Proteins,Number_Layers,Number_Diseases,SupraBipartiteMatrix,lambda){
  
  Transition_Disease_Protein <- Matrix(0,nrow=Number_Diseases,ncol=Number_Proteins*Number_Layers,sparse = TRUE)
  
  colnames(Transition_Disease_Protein) <- rownames(SupraBipartiteMatrix)
  rownames(Transition_Disease_Protein) <- colnames(SupraBipartiteMatrix)
  
  Row_Sum_Bipartite <- rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)
  
  for (i in 1:(Number_Proteins*Number_Layers)){
    if (Row_Sum_Bipartite[i] != 0){
      Transition_Disease_Protein[,i] <- (lambda*SupraBipartiteMatrix[i,])/Row_Sum_Bipartite[i]
    }
  }
  return(Transition_Disease_Protein)
}

# 14.3.-Multiplex intra-transition Matrix. It's very slow... Transform to c++
################################################################################################################
get.transition.multiplex <- function(Number_Proteins,Number_Layers,lambda,SupraAdjacencyMatrix,SupraBipartiteMatrix){
  
  Transition_Multiplex_Network <- Matrix(0,nrow=Number_Proteins*Number_Layers,ncol=Number_Proteins*Number_Layers,sparse = TRUE)
  
  rownames(Transition_Multiplex_Network) <- rownames(SupraAdjacencyMatrix)
  colnames(Transition_Multiplex_Network) <- colnames(SupraAdjacencyMatrix)
  
  Col_Sum_Multiplex <- colSums(SupraAdjacencyMatrix,na.rm=FALSE,dims=1, sparseResult=FALSE)
  Row_Sum_Bipartite <- rowSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE) 
  
  for (j in 1:(Number_Proteins*Number_Layers)){
    if(Row_Sum_Bipartite[j] != 0){
      Transition_Multiplex_Network[,j] <- ((1-lambda)*SupraAdjacencyMatrix[,j]) /Col_Sum_Multiplex[j]
    } else {
      Transition_Multiplex_Network[,j] <- SupraAdjacencyMatrix[,j] /Col_Sum_Multiplex[j]
    }
  }
  return(Transition_Multiplex_Network)
}

# 14.4.-DiseaseSimilarity intra-transition Matrix.
################################################################################################################
get.transition.disease <- function(Number_Diseases,lambda,AdjMatrix,SupraBipartiteMatrix){
  
  Transition_Disease_Network <- Matrix(0,nrow=Number_Diseases,ncol=Number_Diseases,sparse = TRUE)
  
  rownames(Transition_Disease_Network) <- rownames(AdjMatrix)
  colnames(Transition_Disease_Network) <- colnames(AdjMatrix)
  
  Col_Sum_Disease <- colSums (AdjMatrix,na.rm=FALSE,dims=1, sparseResult=FALSE)
  Col_Sum_Bipartite <- colSums (SupraBipartiteMatrix, na.rm = FALSE, dims = 1,sparseResult = FALSE)
  
  for (j in 1:Number_Diseases){
    if(Col_Sum_Bipartite[j] != 0){
      Transition_Disease_Network[,j] <- ((1-lambda)*AdjMatrix[,j]) /Col_Sum_Disease[j]
    } else {
      Transition_Disease_Network[,j] <- AdjMatrix[,j] /Col_Sum_Disease[j]
    }
  }
  return(Transition_Disease_Network)
}
