################################################################################################################
################################################################################################################
####
#### Generate_Pathway_Network: # Protocol to Build Interaction Network from Pathway public Data.
#### 
################################################################################################################
################################################################################################################
####
#### 1. NAME: Generate_Pathway_Network.R 
#### 2. CONTENTS: Protocol to Build Interaction Network from pathway data retrieved from GRAPHITE package: 
####    fetchs the pathways in several databases (reactome, Biocarta, Kegg, Panther and nic) and from that
####    information links proteins. 
#### 3. CREATED: 23/11/2016 by Alberto Valdeolivas.
####
################################################################################################################
################################################################################################################
#### 4. DESCRIPTION:Protocol to Build Interaction Network from pathway data retrieved from GRAPHITE package: 
####    fetchs the pathways in several databases (reactome, Biocarta, Kegg, Panther and nic) and from that
####    information links proteins. 
####    - GRAPHITE Package: Sales G, Calura E, Romualdi C. 2014. 
####    graphite: GRAPH interaction from pathway topological environment. graphite_1.20.1 (Version)
####    - BioMart R package version: biomaRt_2.30.0
################################################################################################################
################################################################################################################
#### 5. OUTPUT: A network file, cotanining two columns with HGNC symbols. Genes in the same record are 
#### linked in the network. The file is named as Pathway_date.gr, where date is replaced for the date 
#### of execution of the script.
################################################################################################################
################################################################################################################
#### 6. EXAMPLE OF EXECUTION: Rscript Generate_Pathway_Network.R 
################################################################################################################
################################################################################################################

rm(list=ls())
################################################################################################################
## 1.- We load the common network generation functions, and we install/load the packages needed. 
################################################################################################################
source("Network_Generation_Utils.R")
CRAN.packages <- c("igraph")
bioconductor.packages <- c("graphite","org.Hs.eg.db")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)

# We create a directoy called networks_files to save the files generated.
data.dir <- "networks_files"
if (!file.exists(data.dir)) {
  dir.create(data.dir)
}

# We get the date of execution. Then we will paste it to the name of the network. 
date <- Sys.Date()

################################################################################################################
## 2.- We recover the different pathways DB. 
################################################################################################################
pathwayDatabases() 

######## 2.1.-KEGG
# We create Kegg interaction Network using the GetPathwayData funciton. 
################################################################################################################
KeggRaw <- GetPathwayData ("kegg","hsapiens")

# Build network.
Kegg <- build.network(KeggRaw)

# ecount(Kegg)
# vcount(Kegg)
# graph_name <- paste("networks_files/Kegg_",date,".gr",sep="",collapse = "")
# write.graph(Kegg, graph_name, format="ncol", weights=NULL)

######## 2.2.-BIOCARTA
# We create Biocarta Network using the GetPathwayData funciton. 
################################################################################################################
BiocartaRaw <- GetPathwayData ("biocarta","hsapiens")

# Build network.
Biocarta <- build.network(BiocartaRaw)

# ecount(Biocarta)
# vcount(Biocarta)
# graph_name <- paste("networks_files/Biocarta_",date,".gr",sep="",collapse = "")
# write.graph(Biocarta, graph_name, format="ncol", weights=NULL)

######## 2.3.-PID Pathway Interaction Database (NCI)
# We create PID Network using the GetPathwayData funciton. 
################################################################################################################
nciRaw <- GetPathwayData ("nci","hsapiens")

# Check Network
nci <- build.network(nciRaw)

# ecount(nci)
# vcount(nci)
# graph_name <- paste("networks_files/nci_",date,".gr",sep="",collapse = "")
# write.graph(nci, graph_name, format="ncol", weights=NULL)

######## 2.4.-Reactome Pathway Interaction Database 
# We create Reactome Network using the GetPathwayData funciton. 
################################################################################################################
reactomeRaw <- GetPathwayData ("reactome","hsapiens")

# Check network   
reactome <- build.network(reactomeRaw)

# ecount(reactome)
# vcount(reactome)
# graph_name <- paste("networks_files/reactome_",date,".gr",sep="",collapse = "")
# write.graph(reactome, graph_name, format="ncol", weights=NULL)

######## 2.5.- Panther Pathway Interaction Database 
# We create Panther Network using the GetPathwayData funciton. 
################################################################################################################
######## Panther
PantherRaw <- GetPathwayData ("panther","hsapiens")

# Check network   
Panther<- build.network(PantherRaw)

# ecount(Panther)
# vcount(Panther)
# graph_name <- paste("networks_files/Panther_",date,".gr",sep="",collapse = "")
# write.graph(Panther, graph_name, format="ncol", weights=NULL)

################################################################################################################
## 3.- Fusion of all the pathway databases into one Pathways network and network saving.
################################################################################################################
PathwayRaw <- graph.union(reactome, Biocarta, Kegg, nci,Panther)

Pathway <- igraph::simplify(PathwayRaw)

print(paste("Number of Edges: ", ecount(Pathway),sep="",collapse = NULL))
print(paste("Number of Nodes: ", vcount(Pathway),sep="", collapse = NULL))
graph_name <- paste("networks_files/AllPathways_",date,".gr",sep="",collapse = "")
write.graph(Pathway, graph_name, format="ncol", weights=NULL)

