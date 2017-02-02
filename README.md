# RWR-M and RWR-MH

There are three different elements in this repository: 

1. **/Scripts_and_Files Folder**: It contains the scripts and the files used to obtain the results shown in the article: "Random Walk with Restart on Multiplex and Heterogeneous Biological Networks". 

2. **/RWR-M.zip**: It contains all the files and scripts needed to perform a random walk with restart on a multiplex network (RWR-M), such as described in the article: "Random Walk with Restart on Multiplex and Heterogeneous Biological Networks".

3. **/RWR-MH.zip:** It contains all the files and scripts needed to perform a random walk with restart on a multiplex and heterogeneous network (RWR-MH), such as described in the article: "Random Walk with Restart on Multiplex and Heterogeneous Biological Networks".

## RWR-M 

Random walk with restart on a biological multiplex network integrated by 3 layers accounting for different types of links among genes (PPI, pathways and co-expression). The user provides an initial gene or set of genes, seed nodes, and the algorithm will measure the proximity of every gene within the network to the input seeds. 

#### Installation

All the scripts and files needed to run the algorithm can be obtained by download the compressed file [/RWR-M.zip] (https://github.com/alberto-valdeolivas/RWR-MH/blob/master/RWR-M.zip). Extract the file in the desired location within your computer. The scripts are written in R with some functions developped in C++. Therefore, in order to execute it you need to install [R] (https://cran.r-project.org/doc/manuals/r-release/R-admin.html). 

#### Usage

Once the compressed file is extracted, you should move to this folder. The main script is called RWR-M.R and if you open it you can find a detailed explanation of the process along with instructions about how to execute it, the input required files and the output files that will be created. A case example with all the required files is provided with the downloaded compressed file. To run the algorithm, you need to type in the command line the following order: 

`Rscript RWR-M.R Input_Files/Seeds_Example.txt Input_Files/Parameters_Example.txt Results_Example`

1. **Input_Files/Seeds_Example.txt**: first argument makes reference to the seed genes. It should be a plain text file containing the seeds gene names, (One gene per line and at least one gene), **HGCN gene official symbol required**. In the example provided, TP53 and LMNB1 genes are took as seeds. To use your own seeds you can either modify this file or reference this argument to another plain text file containing your gene names. 

2. **Input_Files/Parameters_Example.txt**: second entry refers to the different parameters of the RWR-M. It should be a plain text file containing the parameters as specified in the example file provided. We strongly recommend to use this file and just modify the numerical values if you want to tune the parameters. The original example file accounts for their default values used in our article. These default values have shown to be very suitable for the method.

3. **Results_Example**: third argument indicates the name you want to provide to the output files. The output files will be generated in the Output_Files folder with the name provided by the user in this argument. Two files with the same name but different extension will be created. The first one is a .txt file containing all the genes of the network and their associated score. The file is sorted by score in such a way that the closer genes to the seeds are located in the top positions. The second file is a .gr network file. It accounts for all the interactions among the seeds and the top k retrieved genes (k is a parameter that we can change in the previous file). This file can be loaded for instance into [Cytoscape](http://www.cytoscape.org/) providing a suitable view of the RWR-M results. 

## RWR-MH 

Random walk with restart on a biological multiplex and heterogeneous network integrated by the 3-layers multiplex network described in the previous paragraph, a diseases-disease similarity network and bipartite interactions between diseases and their associated genes. The user provides an initial gene/disease or set of genes and diseases, seed nodes, and the algorithm will measure the proximity of every gene and disease within the network to the input seeds. 

#### Installation

All the scripts and files needed to run the algorithm can be obtained by download the compressed file [/RWR-MH.zip] (https://github.com/alberto-valdeolivas/RWR-MH/blob/master/RWR-MH.zip). Extract the file in the desired location within your computer. The scripts are written in R with some functions developped in C++. Therefore, in order to execute it you need to install [R] (https://cran.r-project.org/doc/manuals/r-release/R-admin.html). 

#### Usage

Once the compressed file is extracted, you should move to this folder. The main script is called RWR-MH.R and if you open it you can find a detailed explanation of the process along with instructions about how to execute it, the input required files and the output files that will be created. A case example with all the required files is provided with the downloaded compressed file. In fact, the example corresponds to the situation discussed in the paragraph 3.5 of the above-mentioned article. To run the algorithm, you need to type in the command line the following order: 

`Rscript RWR-MH.R Input_Files/Seeds_Example.txt Input_Files/Parameters_Example.txt Results_Example`

1. **Input_Files/Seeds_Example.txt**: first argument makes reference to the seed genes and/or seed diseases. It should be a plain text file containing the seeds gene names, **HGCN gene official symbol required**, and/or the seed **diseases MIM code** (One seed per line and at least one seed). In the example provided, the Wiedemann-Rautenstrauch syndrome, also known as progeorid neonatal syndrome, (MIM code: 264090) is used as the unique seed node. To use your own seeds you can either modify this file or reference this argument to another plain text file containing your gene names. 

2. **Input_Files/Parameters_Example.txt**: second entry refers to the different parameters of the RWR-MH. It should be a plain text file containing the parameters as specified in the example file provided. We strongly recommend to use this file and just modify the numerical values if you want to tune the parameters. The original example file accounts for their default values used in our article. These default values have shown to be very suitable for the method.

3. **Results_Example**: third argument indicates the name you want to provide to the output files. The output files will be generated in the Output_Files folder with a name based on the one provided by the user in this argument. Three files different files will be generated. The first one is a .txt file containing all the genes of the multiplex network and their associated score. The second one is a .txt file containing all the diseases of the disease-disease similarity network and their related score. Both files are sorted by score in such a way that the closer genes/diseases to the seeds are located in the top positions. The third file is a .gr network file. It accounts for all the interactions among the seeds and the top k retrieved genes (k is a parameter that we can change in the previous file). This file can be loaded for instance into [Cytoscape](http://www.cytoscape.org/) providing a suitable view of the RWR-M results. 
