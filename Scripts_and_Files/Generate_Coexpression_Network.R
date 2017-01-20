################################################################################################################
################################################################################################################
####
#### Generate_Coexpression_Network: Protocol to Build Interaction Network from Co-Expression Data
#### 
################################################################################################################
################################################################################################################
####
#### 1. NAME: Generate_Coexpression_Network.R 
#### 2. CONTENTS: Protocol to Build Interaction Network from Co-Expression Data. 
#### 3. CREATED: 23/11/2016 by Alberto Valdeolivas.
####
################################################################################################################
################################################################################################################
#### 4. DESCRIPTION: A Co-Expression Network is generated from RNA-seq expression data retrieved from Protein
#### Atlas webSite: http://www.proteinatlas.org/about/download. The files are named "rna_celline.csv.zip", 
#### and "rna_celline.csv.zip". They correspond to RNA-seq expression levels in 45 cells and 32 tissues.
#### It contains Ensembl genes, corresponding FPKM values (Fragments Per Kilobase Of Exon Per Million Fragments Mapped) 
#### and HGNC canonical names. The correlation threshold is set to 0.7. All pairs of genes that have a correlation 
#### above this threshold (positive or negative) are considered as linked in the Co-Expression network.
#### - Biomart R package version: biomaRt_2.30.0
################################################################################################################
################################################################################################################
#### 5. OUTPUT: A network file, cotanining two columns with HGNC symbols. Genes in the same record are 
#### linked in the network. The file is named as Co-expresion_date.gr, where date is replaced for the date 
#### of execution of the script (Date of generation of the Network and download ofdata from Protein Atlas.)
################################################################################################################
################################################################################################################
#### 6. EXAMPLE OF EXECUTION: Rscript Generate_Coexpression_Network.R 
################################################################################################################
################################################################################################################

rm(list=ls())
################################################################################################################
## 1.- We load the common network generation functions, and we install/load the packages needed. 
################################################################################################################
source("Network_Generation_Utils.R")
CRAN.packages <- c("igraph", "reshape")
bioconductor.packages <- c("biomaRt")
install.packages.if.necessary(CRAN.packages,bioconductor.packages)

# We create a directoy called networks_files to save the files generated.
data.dir <- "networks_files"
if (!file.exists(data.dir)) {
  dir.create(data.dir)
}

# We get the date of execution. Then we will paste it to the name of the network. 
date <- Sys.Date()


################################################################################################################
## 2.- We Fetch RNA-expression data from Protein Atlas. We save Ensemble ID, Gene Name, Sample and FPKM value 
################################################################################################################
## Tissue. 
tissue.url <- "http://www.proteinatlas.org/download/rna_tissue.csv.zip"
data.tissue <- get.data.from.zip(download.if.necessary(tissue.url, data.dir))

## Cell Lines.
celline.url <- "http://www.proteinatlas.org/download/rna_celline.csv.zip"
data.celline <- get.data.from.zip(download.if.necessary(celline.url, data.dir))

## We merge records from tissue and cell lines.
data_Coexpr <- rbind(data.tissue,data.celline)

# Select column of interest (Ensemble ID, Gene.name, Sample, Value)
data2_Coexpr <-data_Coexpr[,c(1:4)]

################################################################################################################
## 3.- We Check that the HGNC symbols of our Co-expression data are in Biomart. (There are some non-coding RNA) 
## We use a Biomart query.
################################################################################################################
ensemble_HGNC <- check.HGNC.symbols(unique(data2_Coexpr$Gene.name))

data3_Coexpr <- data2_Coexpr[which(data2_Coexpr$Gene.name %in% ensemble_HGNC),]

################################################################################################################
## 4.- We transform our data into a matrix in order to calculate the Correlation.
################################################################################################################
## Output a matrix with Gene.name as first column and, in every other column compute the mean Value for each given Sample  
## Then transform "Gene.name"" column into row names
matrix_Coexpr<-cast(data3_Coexpr, Gene.name ~ Sample,mean,value="Value") 
rownames(matrix_Coexpr)<-matrix_Coexpr[,1]
matrix_Coexpr[,1]<-NULL
matrix_Coexpr<-as.matrix(matrix_Coexpr)

### We Compute Spearman correlation (correlation between ranks of values)  
### We need to take the transposed matrix (rows <-> columns)  
### Values are rounded to the 2nd decimal
CR<-round(cor(t(matrix_Coexpr), method="spearman"),2)

################################################################################################################
## 5.- We Select correlation threshold and build a network.
################################################################################################################
### The correlation threshold is set to 0.7. All pairs of genes that have a correlation above this threshold 
### (positive or negative) are considered as linked in the Co-Expression network.  

## Matrix indeces satisfying the criteria.
inds <- which(abs(CR[1:nrow(CR),]) >= 0.7, arr.ind=TRUE)

#retrieve names matching the index
rnames <- rownames(CR)[inds[,1]]
cnames <- colnames(CR)[inds[,2]]

#construct a network from indexes
net<-apply(cbind(rnames, cnames),1, function(x) unname(unlist(x)))
result<-t(net)

# transform matrix into data.frame
Network_Coexpr<-as.data.frame(result,stringsAsFactors = FALSE)

# change column names
colnames(Network_Coexpr) <- c("Gene1", "Gene2")

#transform into matrix
ExprNet<-as.matrix(Network_Coexpr)

################################################################################################################
## 6.- We transform the matrix into a network, remove selfloops and isolated components.
################################################################################################################

### We call the funciton build.network for sanity checks and graph conversion
Expr <- build.network(ExprNet)
print(paste("Number of Edges: ", ecount(Expr),sep="",collapse = NULL))
print(paste("Number of Nodes: ", vcount(Expr),sep="", collapse = NULL))
graph_name <- paste("networks_files/Co-Expression_",date,".gr",sep="",collapse = "")

write.graph(Expr, graph_name, format=c("ncol"),  weights=NULL)
