################################################################################################################
################################################################################################################
####
#### CDF_Plots_Generation_Parameters: # Script that generates a file containing a CDF plot. CDF: Cumulative Distributive
####                        Function. It is the figure associated to the Results 3.4 section. Figure 5
####                                          
################################################################################################################
################################################################################################################
####
#### 1. NAME: CDF_Plots_Generation_Parameters.R 
#### 2. CONTENTS: Protocol to generate the CDF plots of the different options treated on the section Results 3.5
#### 3. CREATED: 13/12/2016 by Alberto Valdeolivas.
####
################################################################################################################
################################################################################################################
#### 4. DESCRIPTION: The script read the results obtained in LOOCV script, where ability to rank genes of RWR
####                 is tested in for different parameters values of RWR-MH. With this data we generate a CDF and
####                 we display it focused in the first ranks (The method is focused on prioritization. Therefore
####                 that is the interesting part of the plot)
####       
################################################################################################################
################################################################################################################
#### 5. OUTPUT: A plot showing the different CDFs for the different multiplex and monoplex approaches.
################################################################################################################
################################################################################################################
#### 6. EXAMPLE OF EXECUTION: Just Run it. Rscript CDF_Plots_Generation_Parameters.R
####                          
################################################################################################################
################################################################################################################
rm(list=ls())

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

################################################################################################################
## 1.- We load the common network generation functions, and we install/load the packages needed. 
################################################################################################################
source("LOOCV_Utils.R")
CRAN.packages <- c("reshape2","plyr","ggplot2")
bioconductor.packages <- c()
install.packages.if.necessary(CRAN.packages, bioconductor.packages)

# We create a directoy called networks_files to save the files generated.
data.dir <- "Results_34"
if (!file.exists(data.dir)) {
  dir.create(data.dir)
}

################################################################################################################
## 2.- We read the ranking results generated by the RWR algorithm in the different situations.
################################################################################################################

### Restart Probability.
Input_Files_restart_01 <- read.table("Results_34/restart/restart_01.txt", header = TRUE, stringsAsFactors = FALSE)
Input_Files_restart_03 <- read.table("Results_34/restart/restart_03.txt", header = TRUE, stringsAsFactors = FALSE)
Input_Files_restart_05 <- read.table("Results_34/restart/restart_05.txt", header = TRUE, stringsAsFactors = FALSE)
Input_Files_restart_07 <- read.table("Results_34/restart/restart_07.txt", header = TRUE, stringsAsFactors = FALSE)
Input_Files_restart_09 <- read.table("Results_34/restart/restart_09.txt", header = TRUE, stringsAsFactors = FALSE)

### delta.
Input_Files_delta_05 <- read.table("Results_34/delta/delta_05.txt", header = TRUE, stringsAsFactors = FALSE)
Input_Files_delta_01 <- read.table("Results_34/delta/delta_01.txt", header = TRUE, stringsAsFactors = FALSE)
Input_Files_delta_09 <- read.table("Results_34/delta/delta_09.txt", header = TRUE, stringsAsFactors = FALSE)

### Lambda.
Input_Files_lambda_05 <- read.table("Results_34/lambda/lambda_05.txt", header = TRUE, stringsAsFactors = FALSE)
Input_Files_lambda_01 <- read.table("Results_34/lambda/lambda_01.txt", header = TRUE, stringsAsFactors = FALSE)
Input_Files_lambda_09 <- read.table("Results_34/lambda/lambda_09.txt", header = TRUE, stringsAsFactors = FALSE)

### eta.
Input_Files_eta_05 <- read.table("Results_34/eta/eta_05.txt", header = TRUE, stringsAsFactors = FALSE)
Input_Files_eta_01 <- read.table("Results_34/eta/eta_01.txt", header = TRUE, stringsAsFactors = FALSE)
Input_Files_eta_09 <- read.table("Results_34/eta/eta_09.txt", header = TRUE, stringsAsFactors = FALSE)

### tau.
Input_Files_tau_1_1_1 <- read.table("Results_34/tau/tau_1_1_1.txt", header = TRUE, stringsAsFactors = FALSE)
Input_Files_tau_1_01_19 <- read.table("Results_34/tau/tau_1_01_19.txt", header = TRUE, stringsAsFactors = FALSE)
Input_Files_tau_1_19_01 <- read.table("Results_34/tau/tau_1_19_01.txt", header = TRUE, stringsAsFactors = FALSE)

################################################################################################################
## 3.- We generate a data frame with all the rankings we are going to take in account. We calculate CDF
################################################################################################################

### Restart Probability.
restart_01 <- Input_Files_restart_01$Global_Ranking
restart_03 <- Input_Files_restart_03$Global_Ranking
restart_05 <- Input_Files_restart_05$Global_Ranking
restart_07 <- Input_Files_restart_07$Global_Ranking
restart_09 <- Input_Files_restart_09$Global_Ranking

### delta.
delta_01 <- Input_Files_delta_01$Global_Ranking
delta_05 <- Input_Files_delta_05$Global_Ranking
delta_09 <- Input_Files_delta_09$Global_Ranking


### lambda.
lambda_01 <- Input_Files_lambda_01$Global_Ranking
lambda_05 <- Input_Files_lambda_05$Global_Ranking
lambda_09 <- Input_Files_lambda_09$Global_Ranking

### eta.
eta_01 <- Input_Files_eta_01$Global_Ranking
eta_05 <- Input_Files_eta_05$Global_Ranking
eta_09 <- Input_Files_eta_09$Global_Ranking

### tau.
tau_1_01_19 <- Input_Files_tau_1_01_19$Global_Ranking
tau_1_1_1 <- Input_Files_tau_1_1_1$Global_Ranking
tau_1_19_01 <- Input_Files_tau_1_19_01$Global_Ranking

################################################################################################################
## 4.- We generate a data frame with all the rankings we are going to take in account. We calculate CDF
## for the different parameters.
################################################################################################################

### Restart Probability.
################################################################################################################
ggdata_restart <- data.frame(restart_01, restart_03,restart_05, restart_07,restart_09)

ggdata_restart <- melt(ggdata_restart)
ggdata_restart <- ddply(ggdata_restart, .(variable), transform, ecd=ecdf(value)(value)) 


cdf_restart <- ggplot(ggdata_restart, aes(x=value, color=variable, linetype=variable)) + 
  stat_ecdf(size=1.75) +
  labs(x="Rank", y="Cumulative Distribution") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","#336600","#009999"),
                     labels=c(expression(r ~ '= 0.1'),expression(r ~ '= 0.3'),expression(r ~ '= 0.5'), expression(r ~ '= 0.7'), 
                              expression(r ~ '= 0.9'))) +
  scale_linetype_manual(guide=FALSE, values=c(3,2,5,6,7),
                        labels=c(expression(r ~ '= 0.1'),expression(r ~ '= 0.3'),expression(r ~ '= 0.5'), expression(r ~ '= 0.7'), 
                                 expression(r ~ '= 0.9'))) +
  guides(color = guide_legend(override.aes = list(color = c("#999999", "#E69F00", "#56B4E9","#336600","#009999"),
                                                  linetype = c(3,2,5,1,4), size= 1.25)))  +
  theme_bw() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"), 
        legend.title=element_blank(),
        legend.text = element_text(size = 14, face = "bold"), 
        legend.position=c(0.8,0.45),
        legend.key.size = unit(1.25, "cm"),
        legend.text.align = 0,
        legend.background = element_rect(colour='black',size=0.75),
        panel.grid.major = element_line(size = 0.5, colour = 'black',linetype='dotted'),  
        panel.border=element_rect(size=0.75,colour='black'))  

cdf_restart <- cdf_restart + coord_cartesian(xlim = c(1, 60), ylim = c(0,0.55))  + annotate("text", x = 5, y = 0.5, label = "A", size = 15)


### delta.
################################################################################################################
ggdata_delta <- data.frame(delta_01, delta_05,delta_09)

ggdata_delta <- melt(ggdata_delta)
ggdata_delta <- ddply(ggdata_delta, .(variable), transform, ecd=ecdf(value)(value)) 
  

cdf_delta <- ggplot(ggdata_delta, aes(x=value, color=variable, linetype=variable)) + 
  stat_ecdf(size=1.75) +
  labs(x="Rank", y="Cumulative Distribution") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"),
                     labels=c(expression(delta ~ '= 0.1'),expression(delta ~ '= 0.5'),expression(delta ~ '= 0.9'))) +
  scale_linetype_manual(guide=FALSE, values=c(3,2,5),
                        labels=c(expression(delta ~ '= 0.1'),expression(delta ~ '= 0.5'),expression(delta ~ '= 0.9')))+
  guides(color = guide_legend(override.aes = list(color = c("#999999", "#E69F00", "#56B4E9"),
                                                  linetype = c(3,2,5), size= 1.25))) + 
  theme_bw() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"), 
        legend.title=element_blank(),
        legend.text = element_text(size = 16, face = "bold"), 
        legend.position=c(0.8,0.4),
        legend.key.size = unit(1.75, "cm"),
        legend.text.align = 0,
        legend.background = element_rect(colour='black',size=0.75),
        panel.grid.major = element_line(size = 0.5, colour = 'black',linetype='dotted'),  
        panel.border=element_rect(size=0.75,colour='black'))  
        
cdf_delta <- cdf_delta + coord_cartesian(xlim = c(1, 60), ylim = c(0,0.55)) + annotate("text", x = 5, y = 0.5, label = "B", size = 15)

                                
################################################################################################################
### lambda.
################################################################################################################
ggdata_lambda <- data.frame(lambda_01, lambda_05,lambda_09)

ggdata_lambda <- melt(ggdata_lambda)
ggdata_lambda <- ddply(ggdata_lambda, .(variable), transform, ecd=ecdf(value)(value)) 


cdf_lambda <- ggplot(ggdata_lambda, aes(x=value, color=variable, linetype=variable)) + 
  stat_ecdf(size=1.75) +
  labs(x="Rank", y="Cumulative Distribution") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"),
                     labels=c(expression(lambda ~ '= 0.1'),expression(lambda ~ '= 0.5'),expression(lambda ~ '= 0.9'))) +
  scale_linetype_manual(guide=FALSE, values=c(3,2,5),
                        labels=c(expression(lambda ~ '= 0.1'),expression(lambda ~ '= 0.5'),expression(lambda ~ '= 0.9')))+
  guides(color = guide_legend(override.aes = list(color = c("#999999", "#E69F00", "#56B4E9"),
                                                  linetype = c(3,2,5), size= 1.25))) + 
  theme_bw() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"), 
        legend.title=element_blank(),
        legend.text = element_text(size = 16, face = "bold"), 
        legend.position=c(0.8,0.4),
        legend.key.size = unit(1.75, "cm"),
        legend.text.align = 0,
        legend.background = element_rect(colour='black',size=0.75),
        panel.grid.major = element_line(size = 0.5, colour = 'black',linetype='dotted'),  
        panel.border=element_rect(size=0.75,colour='black'))  

cdf_lambda<- cdf_lambda + coord_cartesian(xlim = c(1, 60), ylim = c(0,0.55)) + annotate("text", x = 5, y = 0.5, label = "D", size = 15)


################################################################################################################
### eta.
################################################################################################################
ggdata_eta <- data.frame(eta_01, eta_05,eta_09)

ggdata_eta <- melt(ggdata_eta)
ggdata_eta <- ddply(ggdata_eta, .(variable), transform, ecd=ecdf(value)(value)) 


cdf_eta <- ggplot(ggdata_eta, aes(x=value, color=variable, linetype=variable)) + 
  stat_ecdf(size=1.75) +
  labs(x="Rank", y="Cumulative Distribution") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"),
                     labels=c(expression(eta ~ '= 0.1'),expression(eta ~ '= 0.5'),expression(eta ~ '= 0.9'))) +
  scale_linetype_manual(guide=FALSE, values=c(3,2,5),
                        labels=c(expression(eta ~ '= 0.1'),expression(eta ~ '= 0.5'),expression(eta ~ '= 0.9')))+
  guides(color = guide_legend(override.aes = list(color = c("#999999", "#E69F00", "#56B4E9"),
                                                  linetype = c(3,2,5), size= 1.25))) + 
  theme_bw() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"), 
        legend.title=element_blank(),
        legend.text = element_text(size = 16, face = "bold"), 
        legend.position=c(0.8,0.4),
        legend.key.size = unit(1.75, "cm"),
        legend.text.align = 0,
        legend.background = element_rect(colour='black',size=0.75),
        panel.grid.major = element_line(size = 0.5, colour = 'black',linetype='dotted'),  
        panel.border=element_rect(size=0.75,colour='black'))  

cdf_eta <- cdf_eta + coord_cartesian(xlim = c(1, 60), ylim = c(0,0.55)) + annotate("text", x = 5, y = 0.5, label = "E", size = 15)

################################################################################################################
### tau.
################################################################################################################
ggdata_tau <- data.frame(tau_1_01_19, tau_1_1_1,tau_1_19_01)

ggdata_tau <- melt(ggdata_tau)
ggdata_tau <- ddply(ggdata_tau, .(variable), transform, ecd=ecdf(value)(value)) 


cdf_tau <- ggplot(ggdata_tau, aes(x=value, color=variable, linetype=variable)) + 
  stat_ecdf(size=1.75) +
  labs(x="Rank", y="Cumulative Distribution") +
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"),
                     labels=c(expression(tau ~ '= (1, 0.1, 1.9)/L'),expression(tau ~ '= (1, 1, 1)/L'),expression(tau ~ '= (1, 1.9, 0.1)/L'))) +
  scale_linetype_manual(guide=FALSE, values=c(3,2,5),
                     labels=c(expression(tau ~ '= (1, 0.1, 1.9)/L'),expression(tau ~ '= (1, 1, 1)/L'),expression(tau ~ '= (1, 1.9, 0.1)/L'))) +
  guides(color = guide_legend(override.aes = list(color = c("#999999", "#E69F00", "#56B4E9"),
                                                  linetype = c(3,2,5), size= 1.25))) + 
  theme_bw() +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=18,face="bold"), 
        legend.title=element_blank(),
        legend.text = element_text(size = 16, face = "bold"), 
        legend.position=c(0.7,0.4),
        legend.key.size = unit(1.75, "cm"),
        legend.background = element_rect(colour='black',size=0.75),
        legend.text.align = 0,
        panel.grid.major = element_line(size = 0.5, colour = 'black',linetype='dotted'),  
        panel.border=element_rect(size=0.75,colour='black'))  

cdf_tau <- cdf_tau + coord_cartesian(xlim = c(1, 60), ylim = c(0,0.55)) + annotate("text", x = 5, y = 0.5, label = "C", size = 15)



figure_5 <- multiplot(cdf_delta, cdf_lambda, cdf_tau, cdf_eta, cols=2)


figure_5_Final<- multiplot(cdf_restart,figure_5,cols=1)
# dev.copy(png,'../Figures/figure_5.png')
# dev.off()
# ggsave("../Figures/Figure_5.jpg")

grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 2)))
print(cdf_restart, vp = viewport(layout.pos.row =1, layout.pos.col=1))
print(cdf_delta, vp = viewport(layout.pos.row= 2, layout.pos.col=1))
print(cdf_lambda, vp = viewport(layout.pos.row=2, layout.pos.col=2))
print(cdf_tau, vp = viewport(layout.pos.row=3, layout.pos.col=1))
print(cdf_eta, vp = viewport(layout.pos.row=3, layout.pos.col=2))
dev.copy(png,'../Figures_V2/figure_5_V2.png')
dev.off()
# ggsave("../Figures/Figure_5.jpg")
