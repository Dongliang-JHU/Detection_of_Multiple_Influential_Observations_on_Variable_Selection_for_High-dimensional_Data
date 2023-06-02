# Author: Dongliang Zhang (^)
#
# Supervisors: Professor Masoud Asgharian (*)
#              Professor Martin Lindquist (^)
#              Professor Mei-Cheng Wang (^)
#
# Current Affiliation: (*) Department of Mathematics and Statistics, McGill University, Montreal, Quebec, Canada 
#                      (^) Department of Biostatistics, Johns Hopkins University, Baltimore, Maryland, USA  
#
# Title of Project: Detection of Multiple Influential Observations on Variable Selection 
#
# R Script Purpose: Real Operation Platform I (Real Dataset I) 
# 
# Created on  : March 27, 2023 
#
# Modified on : March 27, 2023 
#   
################################################################################ 
################################ Paths ######################################### 
################################################################################ 

# data
dat.path <- "/Users/dongliangzhang/Library/CloudStorage/Dropbox/Research/Detection of Influential Points on Variable Selection/Dongliang Zhang_Influential Diagnostics/Code/Data/"

# support 
support.file.path <- "/Users/dongliangzhang/Library/CloudStorage/Dropbox/Research/Detection of Influential Points on Variable Selection/Dongliang Zhang_Influential Diagnostics/Code/Support/"

#output 
output.file.path <- "/Users/dongliangzhang/Library/CloudStorage/Dropbox/Research/Detection of Influential Points on Variable Selection/Dongliang Zhang_Influential Diagnostics/Code/Output/"

output.dir <- "/Users/dongliangzhang/Library/CloudStorage/Dropbox/Research/Detection of Influential Points on Variable Selection/Dongliang Zhang_Influential Diagnostics/Code/Output/"

#working directory
wd.path <- "/Users/dongliangzhang/Library/CloudStorage/Dropbox/Research/Detection of Influential Points on Variable Selection/Dongliang Zhang_Influential Diagnostics/Code/Development"

setwd(wd.path)

################################################################################ 
############################## Source Files #################################### 
################################################################################ 

source(paste0(support.file.path,"dataGenerate.r"))
source(paste0(support.file.path,"simulation_platform.r"))
source(paste0(support.file.path,"step1_multDetection.r"))
source(paste0(support.file.path,"step1_multCluster.r"))
source(paste0(support.file.path,"step2_filtering.r"))
source(paste0(support.file.path,"normalDecision.r"))
source(paste0(support.file.path,"enet_fit.r"))
source(paste0(support.file.path,"lasso_fit.r"))
source(paste0(support.file.path,"cluster.r"))
source(paste0(support.file.path,"df_measure.r"))
source(paste0(support.file.path,"testsEval.r"))
source(paste0(support.file.path,"binaryClassification.r"))
source(paste0(support.file.path, "real_platform.r"))

################################################################################
############################### Libraries ######################################
################################################################################

library(MASS)
library(mvtnorm)
library(glmnet)
library(ncvreg)
library(scalreg)
library(rje) #is.subset 
library(rlist)
library(cluster)
library(tictoc) 
library(ggplot2)
library(ggpubr)
library(MIP)
library(e1071)
library(SIS)
library(mpath)
library(doParallel)
library(foreach)
library(pROC)
library(ensr)
library(ggpubr)
library(emojifont)

################################################################################ 
#################################### Data ###################################### 
################################################################################ 

#X <- read.csv(paste0(dat.path,"DongliangFeatures.csv"), head = F, sep="")[,-c(471,472,477)]
#Y <- read.csv(paste0(dat.path,"DongliangOutcomes.csv"), head = F, sep="") 

#Xy.dat <- as.data.frame(cbind(Y,X))
#names(Xy.dat)<-c("Y",paste("X",1:ncol(X),sep=""))
#write.csv(Xy.dat, paste0(dat.path,"brainPain.csv"))

dat <- read.csv(paste0(dat.path,"brainPain.csv"))[,-1]

y <- as.vector(dat$Y)

X <- dat[,-1]

n <- nrow(X)
p <- ncol(X)

################################################################################
# identified indices of influential observations 

ClusMIP_LASSO <- c(3, 7, 13, 14, 15, 19, 22, 25, 26, 43, 44, 49,50, 52, 55, 56, 
                   57, 58, 62, 67, 75, 79, 80,81, 103, 115, 122, 123, 124, 127, 
                   145, 147, 157, 158, 159, 164, 165, 169, 170, 172,175, 176, 194)

ClusMIP_enet <- c(31, 49, 50, 52, 55, 58, 76, 79, 115, 122, 124, 158, 169, 170, 172, 194)

ClusMIP_SCAD <- c(31, 43, 49, 52, 55, 57, 62, 76, 80, 86, 103, 121, 133, 134, 
                  145, 146, 157, 159, 164, 165, 173, 177, 181, 182, 193)

ClusMIP_MCP <- c(13, 19, 26, 49, 52, 55, 57, 58, 73, 76, 86, 103, 121, 124, 134, 
                 145, 147, 151, 163, 164, 165, 182)

DF <- c(2, 57, 84, 119, 153, 159, 168, 174, 198)

################################################################################
# Data Processing 

# Assessment I: Prediction Power  

cor.mat <- as.matrix(read.csv("/Users/dongliangzhang/Library/CloudStorage/Dropbox/Research/Detection of Influential Points on Variable Selection/Dongliang Zhang_Influential Diagnostics/Code/Final Result/Real Dataset I/Prediction/cor_mat.csv"))[,-1]

Method <- rep(c("ClusMIP(LASSO)", "ClusMIP(ENET)", "ClusMIP(SCAD)", "ClusMIP(MCP)"), each=3)
Operation <- rep(c("Full", "Removal1", "Removal2"), 4)
Correlation <- as.numeric(t(cor.mat))
df.prediction <- data.frame(Method=Method, Operation=Operation, Correlation=Correlation)
df.prediction$Method <- factor(df.prediction$Method, levels = c("ClusMIP(LASSO)", "ClusMIP(ENET)", "ClusMIP(SCAD)", "ClusMIP(MCP)"))
df.prediction$Operation <- factor(df.prediction$Operation, levels = c("Full", "Removal1", "Removal2"))

# Assessment II: Pain Ratings 

Temperature <- rep(c("44.3","45.3","46.3","47.3","48.3","49.3"), 33)
Label <- rep("Clean", 198)
Label[ClusMIP_SCAD] <- "Influential"
Y <- y 

df.ratings <- data.frame(Ratings=Y, Temperature=Temperature, Label=Label)
df.ratings$Temperature <- factor(df.ratings$Temperature, levels=c("44.3","45.3","46.3","47.3","48.3","49.3"))
df.ratings$Label <- factor(df.ratings$Label, levels=c("Clean","Influential"))

# Assessment III: Multi-dimensional Scaling 

dist.euclidean <- dist(X, method = "euclidean")
scaled <- cmdscale(dist.euclidean, k=2)

v1 <- scaled[,1]
v2 <- scaled[,2]
Label <- rep("Clean", 198)
Label[ClusMIP_MCP] <- "Influential"
df.mdscale <- data.frame(V1=v1, V2=v2, Label=Label)
df.mdscale$Label <- factor(df.mdscale$Label, levels=c("Clean","Influential"))

################################################################################

tmp.size <- 38

#Plot 
jpeg(filename=paste0(output.dir,"RealDataI_Plot_I.jpeg"), width=2400, height=1250) 

txt.ylab <- expression(paste("Correlation between Y and ", hat(Y))) 

txt.xlab.removal0 <- expression('\u2205')
txt.xlab.removal1 <- expression(paste(hat(I)[ClusMIP(Selector)]^'infl')) 
txt.xlab.removal2 <- expression(paste(hat(I)[Common]^'infl')) 

# figure 1 
fig1 <- ggplot(df.prediction, aes(x=Operation, y=Correlation, group=Method)) + 
  geom_line(linetype = "dashed") + 
  geom_point(aes(color=Method), size = 6)  +
  xlab("Data Removal") + ylab(txt.ylab) + 
  #geom_tile() +
  #scale_fill_gradient(low="white", high="red") + 
  theme(axis.title.x = element_text(color="black", size=tmp.size, face="bold", vjust=-0.5),
        axis.title.y = element_text(color="black", size=tmp.size, face="bold"),
        axis.text.x = element_text(color="black",  size=tmp.size, vjust=0),
        axis.text.y = element_text(color="black",  size=tmp.size),
        legend.title = element_text(size=tmp.size,face="bold"),
        #legend.title = element_blank(),  
        legend.text = element_text(size=tmp.size),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = unit(c(3,3,3,3), "lines"),
        legend.position = c(0.7, 0.14))  + 
  scale_y_continuous(breaks = seq(0.5, 1.0, 0.05)) + 
  scale_colour_discrete(name="Detection\nProcedures", labels=c("ClusMIP(LASSO)", 
                                                               "ClusMIP(ENET)", 
                                                               "ClusMIP(SCAD)", 
                                                               "ClusMIP(MCP)")) +
  scale_x_discrete(labels=c(txt.xlab.removal0, txt.xlab.removal1, txt.xlab.removal2))

#figure 2
fig2 <- ggplot(df.ratings, aes(x=Temperature, y=Ratings, group=Label)) + 
  #geom_line(linetype = "dashed") + 
  geom_point(aes(color=Label), size = 6)  +
  xlab(expression("Temperature ("*degree*C*")")) + ylab("Pain Ratings") + 
  #geom_tile() +
  #scale_fill_gradient(low="white", high="red") + 
  theme(axis.title.x = element_text(color="black", size=tmp.size, face="bold", vjust = -0.5),
        axis.title.y = element_text(color="black", size=tmp.size, face="bold"),
        axis.text.x = element_text(color="black",  size=tmp.size, vjust=-0.5),
        axis.text.y = element_text(color="black",  size=tmp.size),
        legend.title = element_text(size=tmp.size,face="bold"),
        #legend.title = element_blank(),  
        legend.text = element_text(size=tmp.size),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = unit(c(3,3,3,3), "lines"),
        legend.position = c(0.19, 0.92)) +
         scale_color_manual(values = c("Clean" = "blue", "Influential" = "red"))

#+ 
  #scale_y_continuous(breaks = seq(0.5, 1.0, 0.05)) + 
  #scale_colour_discrete(name="Detection\nProcedures", labels=c("ClusMIP(LASSO)", 
  #                                                             "ClusMIP(ENET)", 
  #                                                            "ClusMIP(SCAD)", 
  #                                                             "ClusMIP(MCP)")) +
  #scale_x_discrete(labels=c("Full Dataset", txt.xlab.removal1, txt.xlab.removal2))

# figure 3 
fig3 <- ggplot(df.mdscale, aes(x=V1, y=V2, group=Label)) + 
  #geom_line(linetype = "dashed") + 
  geom_point(aes(color=Label), size = 6)  +
  xlab("Dimension 1") + ylab("Dimension 2") + 
  #geom_tile() +
  #scale_fill_gradient(low="white", high="red") + 
  theme(axis.title.x = element_text(color="black", size=tmp.size, face="bold", vjust = -0.5),
        axis.title.y = element_text(color="black", size=tmp.size, face="bold"),
        axis.text.x = element_text(color="black",  size=tmp.size, vjust=-0.5),
        axis.text.y = element_text(color="black",  size=tmp.size),
        legend.title = element_text(size=tmp.size,face="bold"),
        #legend.title = element_blank(),  
        legend.text = element_text(size=tmp.size),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = unit(c(3,3,3,3), "lines"),
        legend.position = c(0.19, 0.92))  +
  scale_color_manual(values = c("Clean" = "blue", "Influential" = "red"))
#+ 
  #scale_y_continuous(breaks = seq(0.5, 1.0, 0.05)) + 
  #scale_colour_discrete(name="Detection\nProcedures", labels=c("ClusMIP(LASSO)", 
  #                                                             "ClusMIP(ENET)", 
  #                                                             "ClusMIP(SCAD)", 
  #                                                             "ClusMIP(MCP)")) +
  #scale_x_discrete(labels=c("Full Dataset", txt.xlab.removal1, txt.xlab.removal2))

final <- ggarrange(fig1, fig2, fig3,  
                   labels = c("(A)", "(B)", "(C)"),
                   ncol = 3, nrow = 1, vjust=+1.5, hjust=-0.60, common.legend = FALSE, 
                   font.label=list(color="black", size=tmp.size, face="bold"))

annotate_figure(final, top = text_grob("", color = "black", face = "bold", 
                                       size = tmp.size))

print(final)

dev.off()



