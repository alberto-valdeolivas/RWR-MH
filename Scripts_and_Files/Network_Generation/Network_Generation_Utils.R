################################################################################################################
################################################################################################################
####
#### Network_Generation_Utils: COMMON FUNCTIONS/APPLICATIONS USED TO GENERATE THE NETWORKS.
#### 
################################################################################################################
################################################################################################################

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
    source("http://bioconductor.org/biocLite.R")
  }
  for (p in bioconductor.packages) {
    if (!require(p, character.only=T)) {
      biocLite(p) 
      library(p, character.only=T)
    }
  }
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 2.- Build a data.frame from a zipped csv file. 	
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
get.data.from.zip <- function (path.to.file, sep=","){		
  internal.file <- tools::file_path_sans_ext(basename(path.to.file))		
  unzipped.file <- unz(path.to.file, internal.file)	
  return(read.table(unzipped.file, header=T, sep=sep))	
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 3.- Download a file from the internet if it does not already exist locally.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#' ### Download a file from the internet if it doesn't already exist locally	
#' If a local file already exists, do not download unless force is TRUE
#' url: string: url where the file is located  	
#' download.dir: string: path to download directory	
#' force: boolean: do we want to override local file with downloaded file	
download.if.necessary <- function (url, download.dir=".", force=F){		
  filename = basename(url)	
  path.to.file = file.path(download.dir, filename)	
  if(!file.exists(path.to.file) | force){ # if file does not exists or force download	
    download.file(url, path.to.file, method="curl", mode='w')	
  }	
  return(path.to.file)	
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 4.- We check if the HGNC symbols that we will use to generate the networks
####  exist in Biomart.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
check.HGNC.symbols <- function(GeneNames){

ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")

ensemble_HGNC <- getBM(attributes=c('hgnc_symbol'), 
                       filters='hgnc_symbol', values=GeneNames, mart=ensembl)

return(ensemble_HGNC$hgnc_symbol)
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 5.- We built a network 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#' Build an igraph network	
#' Take df with 2 columns, format "ncol" (one line per interaction, the 2 interactors are separated by a tab)
#' returns a simplified network as an igraph object.	
build.network <- function (raw.network) {	
  net <- igraph::graph.edgelist(raw.network, directed=F)	
  # simplify
  net <- igraph::simplify(net,remove.multiple = TRUE, remove.loops = TRUE)
  
  # We additionally delete isolated nodes.
  isolates <- which(degree(net, mode = c("all")) == 0) - 1
  net <- delete.vertices(net, names(isolates))
  
  return (net) 	
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# 6.- We built a network 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# Take 2 specific columns, remove lines containing incomplete data
clean.network <- function(net) {
  net <- net[, c("Symbol.A", "Symbol.B")]
  net[net == ""] <- NA
  net[net == "-"] <- NA
  return(na.omit(net))
}


#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### 7.- Function to fetch pathway data from pathway databases using Graphite.  
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
## Uses the graphite package (Sales G, Calura E and Romualdi C (2014). graphite: 
## GRAPH Interaction from pathway Topological Environment. R package version 1.12.0.). 

GetPathwayData <- function (db.name, db.specie){
  db <- pathways(db.specie, db.name)
  OneDb.HGNC <- data.frame()
  
  for (i in 1:length(db)) {
    current.pathway <- (convertIdentifiers(db[[i]],"SYMBOL"))
    current.pathway.info <- try(data.frame(src = current.pathway@edges$src,
                                           dest = current.pathway@edges$dest,
                                           direction = current.pathway@edges$direction,
                                           type = current.pathway@edges$type,
                                           name = current.pathway@title),silent=TRUE) 
    OneDb.HGNC <- rbind(OneDb.HGNC,current.pathway.info)
  }
  
  OneDb.HGNC<-OneDb.HGNC[complete.cases(OneDb.HGNC),] #remove NA rows
  
  #transform the directed edge into undirected edge
  OneDb.HGNC[OneDb.HGNC == "directed"] <- "undirected" 
  
  #create a sorted labelled per edge
  OneDb.HGNC$name1 <- pmin(as.vector(OneDb.HGNC$src), as.vector(OneDb.HGNC$dest))
  OneDb.HGNC$name2 <- pmax(as.vector(OneDb.HGNC$src), as.vector(OneDb.HGNC$dest))
  OneDb.HGNC$inter.label <- paste(OneDb.HGNC$name1,OneDb.HGNC$name2, sep="_")
  
  #select the unique inter.label = unique edge
  OneDb.HGNC <- subset(OneDb.HGNC,!(duplicated(OneDb.HGNC$inter.label)))
  
  # select the 2 columns of interest
  NetworkRaw<-as.matrix(OneDb.HGNC[,c(1,2)])
  colnames(NetworkRaw) <- c("SymbolA", "SymbolB")
  
  return(NetworkRaw)
  
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### 8.- Function to calculate the Information content similarity between 
####     two phenotypes.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
# t- One phenotype.
# s- Other phenotype.
# All_terms_Freq - Information content for all the phenotypes based on their
# frequency in the ontology graph
# HPO_GO <- R Object containing the structure of the ontology. (HPO)
Two_Term_Similarity <- function (t,s,All_terms_IC,Graph_GO ) {
  
  ## For the similarity we have to take in account ancestors.
  ancestors_t <- get.ancestors(Graph_GO,t)
  ancestors_s <- get.ancestors(Graph_GO,s)
  
  intersection_ancestors <- ancestors_t[which(ancestors_t %in% ancestors_s)]
  
  # We return the maximum of the information content for each intersected term found
  return(max(All_terms_IC[which(names(All_terms_IC) %in% intersection_ancestors)]))
} 

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### 9.- Function to calculate the the similarity between each pair of disease.
#### It's adapted two run in Parallelization. 
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
Disease_Similarity_Par <- function(x, list_disease){
  
  DiseaseA <- list_disease[[x[[1]]]]
  DiseaseB <- list_disease[[x[[2]]]]
  length_A <- length(DiseaseA)
  length_B <- length(DiseaseB)
  
  Sum_A <- 0
  
  for (i in 1:length_A){
    sim_value <- numeric(length_A * length_B)
    for (j in 1:length_B){
      sim_value[i*j] <- Two_Term_Similarity(DiseaseA[i],DiseaseB[j],IC_All_terms,HPO_GO)
    }
    Sum_A <- Sum_A + max(sim_value)
  }
  
  Term_A <- Sum_A/(2*length_A)
  
  Sum_B <- 0
  
  for (i in 1:length_B){
    sim_value <- numeric(length_B * length_A)
    for (j in 1:length_A){
      sim_value[j*i] <- Two_Term_Similarity(DiseaseB[i],DiseaseA[j],IC_All_terms,HPO_GO)
    }
    Sum_B <- Sum_B + max(sim_value)
  }
  
  Term_B <- Sum_B/(2*length_B)
  
  return((Term_A + Term_B))
}

#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
#### 10.- Function to generate a KNN disease-disease similarity graph
#### from the data frame containing the similarity score between each pair of
#### diseases.
#### #### #### #### #### #### #### #### #### #### #### #### #### #### #### ####
get.knn.graph.from.similarity <- function(df_similarity,k){
  
  Diseases_unique <- (unique(c(df_similarity[,1], df_similarity[,2])))
  m <- length(Diseases_unique)
  Disease_Sim_Matrix <- matrix(0,nrow=m,ncol=m,dimnames = list(Diseases_unique,Diseases_unique))
  
  for (i in 1:m){
    current_disease <- Diseases_unique[i]
    colA <- df_similarity[which(df_similarity[,1] == current_disease),c(2,3)]
    colB <- df_similarity[which(df_similarity[,2] == current_disease),c(1,3)]
    colnames(colA) <- c("Disease","Score")
    colnames(colB) <- c("Disease","Score")
    Scores_Current_Diseases <- rbind.data.frame(colA, colB)
    Scores_Current_Diseases_k <- head(Scores_Current_Diseases[with(Scores_Current_Diseases, order(-Score)), ],k)
    Disease_Sim_Matrix[i, as.character(Scores_Current_Diseases_k$Disease)]   = Scores_Current_Diseases_k$Score 
    Disease_Sim_Matrix[as.character(Scores_Current_Diseases_k$Disease), i]   = Scores_Current_Diseases_k$Score 
  }
  
  disease_network <- graph_from_adjacency_matrix(Disease_Sim_Matrix, mode = c("undirected"), weighted = TRUE, diag = FALSE,
                                                 add.colnames = NULL)
  return(disease_network)
}
