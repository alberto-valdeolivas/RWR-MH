################################################################################################################
################################################################################################################
####
#### Generate_Disease_Network_1:  # First part of the Protocol to Build a Disease-Disease Similarity Network 
####                               from HPO. In this part we calculate similarity between every pair of diseases.
#### 
################################################################################################################
################################################################################################################
####
#### 1. NAME: Generate_Disease_Network_1.R 
#### 2. CONTENTS: Protocol to Build a Disease-Disease Similarity Network from Human Phenotype Ontology (HPO)  
####    data. We calculate the similarty among each pair of diseases based on the information content of 
####    their phenotype ontology (Generate_Disease_Network_1). 
####    Then, that information will be used to build a KNN similarity graph where 
####    each disease is connected to its five nearest neighbours (Generate_Disease_Network_2).
#### 3. CREATED: 24/11/2016 by Alberto Valdeolivas.
####
################################################################################################################
################################################################################################################
#### 4. DESCRIPTION: Protocol to Build a Disease-Disease Similarity Network from Human Phenotype Ontology (HPO)  
####    data. We calculate the similarty among each pair of diseases based on the information content of 
####    their phenotype ontology. A KNN graph is constructed according to their similarity scores, and each
####    disease is connected to its five nearest neighbours.
####    - HPO Version: phenotype_annotation.tab: Compilar #1233 (13-ene-2016 13:18:26)
####    - HPO Ontolohy Graph: "http://purl.obolibrary.org/obo/hp.obo" (Downloaded on 24/11/2016)
################################################################################################################
################################################################################################################
#### 5. OUTPUT: A file containing the similarity score between every pair of diseases from HPO. That file
####    will be after used to generate the disease-disease similarity network in the script 
####    Generate_Disease_Network_2.
################################################################################################################
################################################################################################################
#### 6. EXAMPLE OF EXECUTION: Rscript Generate_Disease_Network_1.R 
################################################################################################################
################################################################################################################

rm(list=ls())
################################################################################################################
## 1.- We load the common network generation functions, and we install/load the packages needed. 
################################################################################################################
source("Network_Generation_Utils.R")
CRAN.packages <- c("igraph","hpoPlot","parallel")
bioconductor.packages <- c("Rgraphviz")
install.packages.if.necessary(CRAN.packages, bioconductor.packages)

# We create a directoy called networks_files to save the files generated.
data.dir <- "networks_files"
if (!file.exists(data.dir)) {
  dir.create(data.dir)
}

# We get the date of execution. Then we will paste it to the name of the network. 
date <- Sys.Date()

################################################################################################################
## 2.- We read disease-phenotype relations from HPO. 
################################################################################################################
# We have to read phenotype_annotation.tab from HPO. 
phenotype_disease_file <- "http://compbio.charite.de/jenkins/job/hpo.annotations/lastStableBuild/artifact/misc/phenotype_annotation.tab"
HPO_DISEASE_PHENOTYPE <- read.csv(download.if.necessary(phenotype_disease_file, data.dir), header=F, sep="\t", 
                         stringsAsFactors = FALSE)
colnames(HPO_DISEASE_PHENOTYPE) <- c("DB","DB_Object_ID","DB_Name","Qualifier","HPO_ID","DB_Reference",
                                     "Evidence_code","Onset_modifier","Frequency_modifier","With","Aspect",
                                     "Synonym","Date","Assigned_by")

## I am just going to keep data coming from OMIM. Additionally we are going to create a a data frame just with
## the IDs of the disease and its associated phenotypes.

HPO_DISEASE_PHENOTYPE_OMIM <- HPO_DISEASE_PHENOTYPE[HPO_DISEASE_PHENOTYPE$DB == "OMIM",]

HPO_DISEASEID_PHENOID <- data.frame(DISEASE_ID = HPO_DISEASE_PHENOTYPE_OMIM$DB_Object_ID ,
                                    PHENOTYPE_ID = HPO_DISEASE_PHENOTYPE_OMIM$HPO_ID
                                    ,stringsAsFactors = FALSE)

HPO_DISEASEID_PHENOID <- HPO_DISEASEID_PHENOID[with(HPO_DISEASEID_PHENOID, order(DISEASE_ID, PHENOTYPE_ID)), ]

################################################################################################################
## 3.- We read the HPO ontology Graph.
################################################################################################################
## We read the OBO file (ONTOLOGY FOMART) --> http://human-phenotype-ontology.github.io/downloads.html) 
hpo_graph <- "networks_files/hp.obo"
HPO_GO <- get.ontology(hpo_graph,qualifier = "HP")

################################################################################################################
## 4.- Dealing with the desease-phenotype in relation with the Ontology graph. 
################################################################################################################

#### 4.1. We convert our data frame to a list. Keys are the diseases containing their associated phenotypes.
################################################################################################################
list_disease_phenos <- list()
a <- 1
i <- 1

while (i < nrow(HPO_DISEASEID_PHENOID)){
  current_disease <- HPO_DISEASEID_PHENOID$DISEASE_ID[i]
  current_phenotpyes <- HPO_DISEASEID_PHENOID$PHENOTYPE_ID[which(HPO_DISEASEID_PHENOID$DISEASE_ID == current_disease)]
  i <- i + length(current_phenotpyes)
  
  ### For each disease we get the minimal set of HPO terms. (See clean.terms)
  current_phenotpyes <- clean.terms(HPO_GO,current_phenotpyes) 
  
  list_disease_phenos[[a]] <- current_phenotpyes
  names(list_disease_phenos)[a] <- current_disease 
  a <- a + 1
}

#### 4.2.With the frequency o all the terms in the ontology we can calculate their Information Content.
################################################################################################################

All_Terms_Frequencies <- get.term.frequencies(HPO_GO, HPO_DISEASEID_PHENOID$PHENOTYPE_ID)
IC_All_terms <- -log(All_Terms_Frequencies)

################################################################################################################
## 5.- SIMILARITY BETWEEN DISEASES. Steps to get the similarity between two diseases.
################################################################################################################
# We prepare the execution in paralel. 
no_cores <- detectCores() 
cl <- makeCluster(no_cores)
clusterExport(cl, c("Two_Term_Similarity","HPO_GO","IC_All_terms"))
clusterEvalQ(cl, library(hpoPlot))

# We generate a data frame which contains all the pairs of possible diseases.
# Very heavy execution. It takes very long. To run in a powerfull computer. 
df <- as.data.frame(t(combn(names(list_disease_phenos), 2)), stringAsFactors = FALSE)
colnames(df) <- c("DiseaseA","DiseaseB")


Similarity_Vector <- parApply(cl,df,1, Disease_Similarity_Par, list_disease_phenos)
stopCluster(cl)

## We write the results. 
All_Diseases_Sim <- cbind(df,Similarity_Vector)
write.table(All_Diseases_Sim,file="All_Diseases.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)




