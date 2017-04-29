################################################################################################################
################################################################################################################
####
#### Generate_Disease_Network_2:  # Second part of the Protocol to Build a Disease-Disease Similarity Network 
####                               from HPO. In this part we generate the KNN graph from the similarity between
####                               each pair of disaeases generated in Generate_Disease_Network_1.
#### 
################################################################################################################
################################################################################################################
####
#### 1. NAME: Generate_Disease_Network_2.R 
#### 2. CONTENTS: Protocol to Build a Disease-Disease Similarity Network from Human Phenotype Ontology (HPO)  
####    data. We calculate the similarty among each pair of diseases based on the information content of 
####    their phenotype ontology (Generate_Disease_Network_1). 
####    Then, that information will be used to build a similarity graph where 
####    each disease is connected to its five nearest neighbours (Generate_Disease_Network_2).
#### 3. CREATED: 29/11/2016 by Alberto Valdeolivas.
####
################################################################################################################
################################################################################################################
#### 4. DESCRIPTION: Protocol to Build a Disease-Disease Similarity Network from Human Phenotype Ontology (HPO)  
####    data. We calculate the similarty among each pair of diseases based on the information content of 
####    their phenotype ontology. A graph is constructed according to their similarity scores, and each
####    disease is connected to its five nearest neighbours.
####    - HPO Version: phenotype_annotation.tab: Compilar #1233 (13-ene-2016 13:18:26)
####    - HPO Ontolohy Graph: "http://purl.obolibrary.org/obo/hp.obo" (Downloaded on...)
################################################################################################################
################################################################################################################
#### 5. OUTPUT: A network file, cotanining two columns with diseases HPO ID. Diseases in the same record are 
#### linked in the network. The file is named as Disease_Similarity_date.gr, where date is replaced for the date 
#### of execution of the script.
################################################################################################################
################################################################################################################
#### 6. EXAMPLE OF EXECUTION: Rscript Generate_Disease_Network_2.R 
################################################################################################################
################################################################################################################

rm(list=ls())
################################################################################################################
## 1.- We load the common network generation functions, and we install/load the packages needed. 
################################################################################################################
source("Network_Generation_Utils.R")
CRAN.packages <- c("igraph")
bioconductor.packages <- c("PSICQUIC")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)

# We create a directoy called networks_files to save the files generated.
data.dir <- "networks_files"
if (!file.exists(data.dir)) {
  dir.create(data.dir)
}

# We get the date of execution. Then we will paste it to the name of the network. 
date <- Sys.Date()

################################################################################################################
## 2.- We read the file generated in Generate_Disease_Network_1.R and we convert into a KNN disease- 
## disease similarity graph (K = 5)
################################################################################################################
All_Diseases <- read.table("All_Diseases.txt", sep="\t", dec=".", stringsAsFactors = FALSE, header = TRUE)

Disease_Similarity_KNN_5 <- get.knn.graph.from.similarity(All_Diseases,k = 5)

################################################################################################################
## 3.- We simplify the network and we save it. 
################################################################################################################
Disease_Similarity_KNN_5 <- igraph::simplify(Disease_Similarity_KNN_5)

print(paste("Number of Edges: ", ecount(Disease_Similarity_KNN_5),sep="",collapse = NULL))
print(paste("Number of Nodes: ", vcount(Disease_Similarity_KNN_5),sep="", collapse = NULL))
graph_name <- paste("networks_files/DiseaseSimilarity_",date,".gr",sep="",collapse = "")
write.graph(Disease_Similarity_KNN_5, graph_name, format="ncol", weights=NULL)



