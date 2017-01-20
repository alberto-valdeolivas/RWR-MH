################################################################################################################
################################################################################################################
####
#### Generate_PPI_Network: # Protocol to Build Interaction Network from PPI Data.
#### 
################################################################################################################
################################################################################################################
####
#### 1. NAME: Generate_PPI_Network.R 
#### 2. CONTENTS: Protocol to Build Interaction Network from PPI data retrieved from PSICQUIC: fetchs  
####    protein-protein interaction data from public databases. Additional interactions are took from
####    CCSB Interactome Database.
#### 3. CREATED: 23/11/2016 by Alberto Valdeolivas.
####
################################################################################################################
################################################################################################################
#### 4. DESCRIPTION: The goal is to construct a protein-protein interaction network containing binary direct 
### physical protein interactions.  Two sources are used : 
### 1.- * The PSICQUIC portal, which allow to retrieve interaction data from different databases participating 
### to the IMEX consortium. These databases describe the interactions in the PSI-MI format, allowing to select 
### direct interactions between human proteins.  
### Reference: del-Toro, N et al.(2013). A new reference implementation of the PSICQUIC web service. 
### Nucleic Acids Research, 41(Web Server issue), W601–6. R package version: PSICQUIC_1.12.0. biomaRt_2.30.0
### 2.- * The Center for Cancer Systems Biology (CCSB) generates human interactome data. They have released in 
### particular the dataset HI-II-2014 containing validated yeast 2-hybrid interactions between human proteins.  
### Reference: Rolland, T. et al. (2014). A Proteome-Scale Map of the Human Interactome Network. 
### Cell, 159(5), 1212–1226. doi:10.1016/j.cell.2014.10.050. 
###
### We will merge the binary interactions described in each approach. 
### - Biomart R package version biomaRt_2.30.0. 
###
################################################################################################################
################################################################################################################
#### 5. OUTPUT: A network file, cotanining two columns with HGNC symbols. Genes in the same record are 
#### linked in the network. The file is named as PPI_date.gr, where date is replaced for the date 
#### of execution of the script (Date of download from PSICQUIC and CCSB.)
################################################################################################################
################################################################################################################
#### 6. EXAMPLE OF EXECUTION: Rscript Generate_PPI_Network.R 
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
## 2.- We Fetch the PSICQUIC databases to obtain human protein-protein interactions. 
################################################################################################################
psicquic <- PSICQUIC()

# Select databases
DB<-c("BioGrid","InnateDB-IMEx","IntAct","MatrixDB","MBInfo",
      "MINT","Reactome","Reactome-FIs","UniProt")
# We retrieve the human interactions from the DB defined above.
tbl.big <- PSICQUIC::interactions(psicquic, species="9606",provider=DB, type="direct interaction")

# Mapping gene names to HGNC  
# Different gene name formats are used in the different databases. Everything is converted to HGNC gene symbol. 
tbl.big <- addGeneInfo(IDMapper("9606"),tbl.big) 

# Extract columns of interest (protein A / protein B) and construct binary network 
# Replace empty cells by NA  
# Remove incomplete lines
colnames(tbl.big)[names(tbl.big) == "A.name"] <- "Symbol.A"
colnames(tbl.big)[names(tbl.big) == "B.name"] <- "Symbol.B"
NetHGNC <- clean.network(tbl.big)

################################################################################################################
## 3.- We take the HUMAN PPI data from CCSB interactome database 
################################################################################################################
# Download the data file from "http://interactome.dfci.harvard.edu/H_sapiens/index.php" 
# Caution : column headers contain many spaces that need to be removed manually.
interactome.file <- "http://interactome.dfci.harvard.edu/H_sapiens/download/HI-II-14.tsv"
CCSB <- read.table(download.if.necessary(interactome.file, data.dir), header=T, sep="\t", as.is=T)

# We Select column of interest, replace empty cells by NA, remove incomplete lines
CCSB <- clean.network(CCSB)

################################################################################################################
## 4.- We merge Networks coming from both data sources and we write out the resulting network.
################################################################################################################
######### Merge the 2 PPI Networks ##################################
net <-rbind (CCSB, NetHGNC)

### We call the funciton build.network for sanity checks and graph conversion
net <- build.network(as.matrix(net))
print(paste("Number of Edges: ", ecount(net),sep="",collapse = NULL))
print(paste("Number of Nodes: ", vcount(net),sep="", collapse = NULL))
graph_name <- paste("networks_files/PPI_",date,".gr",sep="",collapse = "")

## We save the network Print network  
write.graph(net, graph_name, format="ncol", weights=NULL)
