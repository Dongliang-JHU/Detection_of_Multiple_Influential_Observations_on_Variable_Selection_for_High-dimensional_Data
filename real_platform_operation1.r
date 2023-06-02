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
# Created on  : August 24, 2022 
#
# Modified on : August 24, 2022 
#               August 25, 2022 
#               December 09, 2022 
#               March 16, 2023
#               March 23, 2023 
#   
##########################################################################################
########################################## Paths ######################################### 
##########################################################################################

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

##########################################################################################
###################################### Source Files ###################################### 
##########################################################################################

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

##########################################################################################
######################################## Libraries #######################################
##########################################################################################

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

################################################################################ 
########################################## Data ################################ 
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
# clustering 

Xy.dat <- as.data.frame(cbind(X, y))
names(Xy.dat)<-c(paste("X", 1:p, sep=""),"Y")

clean.pos.screening <- cluster.est(df=Xy.dat, K.vec=2) 

infl.pos.screening <- setdiff(1:nrow(X), clean.pos.screening)

label <- rep("Clean", nrow(X))
label[infl.pos.screening] <- "Influential"

pos.dat <- data.frame(pos=1:nrow(X), Group=label)
#df.time$Methods <- factor(df.time$Methods, levels = method.names)
#df.time$Proportion <- factor(df.time$Proportion, levels = infl.prop)

jpeg(filename=paste0("trial.jpeg"), width=2000, height=1600) 
ggplot(pos.dat, aes(x=Group, y=pos, group=Group)) +
  #geom_line(linetype = "dashed") + 
  #geom_point(aes(shape=Methods, color=Proportion), size = 4)  +
  geom_point(aes(color=Group), size = 4)  +
  xlab("") + ylab("Indices") + 
  theme(axis.title.x = element_text(color="black", size=30, face="bold"),
        axis.title.y = element_text(color="black", size=30, face="bold"),
        axis.text.x = element_text(color="black", size=30),
        axis.text.y = element_text(color="black", size=30),
        legend.title = element_text(size=30),
        legend.text = element_text(size=30),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = unit(c(3,3,3,3), "lines")) 
dev.off()
  #+ 
  #scale_y_continuous(breaks = seq(0, 1.0, 0.1)) + 
  #coord_cartesian(ylim = c(0, 1.0)) +   
  #scale_colour_discrete(name="Proportion", labels=infl.prop) + 
  #scale_shape_discrete(name="Methods", labels=method.names)

length(infl.pos.screening)/198

length(clean.pos.screening)/198


################################################################################ 
# detection 

output <- real.detection(X=as.matrix(X), y=y, cluster.type="K.means", 
                         detection.switch=c(1,1,1,1,1,1,1,0), 
                         alpha=0.05, K.vec=2, n_folds=10, n_subset=30, subset.size=nrow(X)/2, family="gaussian")

output2 <- real.detection(X=as.matrix(X), y=y, cluster.type="K.means", 
                          detection.switch=c(0,0,0,1,1,1,1,0), 
                          alpha=0.05, K.vec=2, n_folds=10, n_subset=30, subset.size=nrow(X)/2, family="gaussian")

# total output 
list.save(output, file = paste0(output.file.path, "output_brainPain.RData"))
#C <- list.load(paste0(output.file.path, "output_", iter[M], ".RData"))

write.csv(output$dec.mat, file = paste0(output.file.path, "dec_brainPain.csv"))
write.csv(output$time.mat, file = paste0(output.file.path, "time_brainPain.csv"))
write.csv(output$infl.pos.cluster, file = paste0(output.file.path, "inflPosCluster_brainPain.csv"))
list.save(output$infl.pos, file = paste0(output.file.path, "inflPos_brainPain.RData"))

ClusMIP_LASSO <- c(3, 7, 13, 14, 15, 19, 22, 25, 26, 43, 44, 49,50, 52, 55, 56, 57, 
58, 62, 67, 75, 79, 80,81, 103, 115, 122, 123, 124, 127, 145, 147,
157, 158, 159, 164, 165, 169, 170, 172,175, 176, 194)


ClusMIP_enet <- c(31, 49, 50, 52, 55, 58, 76, 79, 115, 122, 124, 158, 169, 170, 172, 194)

ClusMIP_SCAD <- c(31, 43, 49, 52, 55, 57, 62, 76, 80, 86, 103,
                  121, 133, 134, 145, 146, 157, 159, 164,
                  165, 173, 177, 181, 182, 193)

ClusMIP_MCP <- c(13, 19, 26, 49, 52, 55, 57, 58, 73, 76, 86,
                 103, 121, 124, 134, 145, 147, 151, 163,
                 164, 165, 182)

DF <- c(2, 57, 84, 119, 153, 159, 168, 174, 198)

Reduce(intersect, list(ClusMIP_LASSO,ClusMIP_enet,ClusMIP_SCAD,ClusMIP_MCP))

################################################################################ 
# refitting 

# LASSO

#full 
lasso.full.original <- cv.glmnet(as.matrix(X), y, alpha=1, family="gaussian", nfolds=10)   

coeff.lasso.full.original <- as.numeric(coef(lasso.full.original, s="lambda.min")[-1]) 

var.lasso.full.original <- which(coeff.lasso.full.original!=0)

#removal part1
y.lasso.removal.part1.original <- y[-ClusMIP_LASSO]
X.lasso.removal.part1.original <- X[-ClusMIP_LASSO,]

lasso.removal.part1.original <- cv.glmnet(as.matrix(X.lasso.removal.part1.original), y.lasso.removal.part1.original, alpha=1, family="gaussian", nfolds=10)

coeff.lasso.removal.part1.original <- as.numeric(coef(lasso.removal.part1.original, s="lambda.min")[-1]) 

var.lasso.removal.part1.original <- which(coeff.lasso.removal.part1.original!=0)

#removal part2  
y.lasso.removal.part2.original <- y[-c(49,52,55)]
X.lasso.removal.part2.original <- X[-c(49,52,55),]

lasso.removal.part2.original <- cv.glmnet(as.matrix(X.lasso.removal.part2.original), y.lasso.removal.part2.original, alpha=1, family="gaussian", nfolds=10)   

coeff.lasso.removal.part2.original <- as.numeric(coef(lasso.removal.part2.original, s="lambda.min")[-1]) 

var.lasso.removal.part2.original <- which(coeff.lasso.removal.part2.original!=0)

################################################################################ 
#scaled LASSO 

#full

slasso.full <- scalreg(as.matrix(X), y, lam0 = NULL, LSE = FALSE) 

coeff.slasso.full <- as.numeric(coef(slasso.full, s="lambda.min")) 

var.slasso.full <- which(coeff.slasso.full!=0)

y.hat.slasso.full <- as.matrix(X) %*% coeff.slasso.full

pred.error.slasso.full <- sum((y-as.vector(y.hat.slasso.full))^2)

cor(y, y.hat.slasso.full)

################################################################################ 
# SCAD 

#full
scad.full.original <- cv.ncvreg(as.matrix(X), y, family="gaussian", penalty="SCAD") 

scad.full.obj.original <- scad.full.original$fit 

coeff.scad.full.original <- as.numeric(scad.full.obj.original$beta[,scad.full.original$min][-1])

var.scad.full.original <- which(coeff.scad.full.original!=0)

#removal part1 
y.scad.removal.part1.original <- y[-ClusMIP_SCAD]
X.scad.removal.part1.original <- X[-ClusMIP_SCAD,]

scad.full.removal.part1.original <- cv.ncvreg(as.matrix(X.scad.removal.part1.original), y.scad.removal.part1.original, family="gaussian", penalty="SCAD") 

scad.full.obj.removal.part1.original <- scad.full.removal.part1.original$fit 

coeff.scad.removal.part1.original <- as.numeric(scad.full.obj.removal.part1.original$beta[,scad.full.removal.part1.original$min][-1])

var.scad.removal.part1.original <- which(coeff.scad.removal.part1.original!=0)

#removal part2 
y.scad.removal.part2.original <- y[-c(49,52,55)]
X.scad.removal.part2.original <- X[-c(49,52,55),]

scad.full.removal.part2.original <- cv.ncvreg(as.matrix(X.scad.removal.part2.original), y.scad.removal.part2.original, family="gaussian", penalty="SCAD") 

scad.full.obj.removal.part2.original <- scad.full.removal.part2.original$fit 

coeff.scad.removal.part2.original <- as.numeric(scad.full.obj.removal.part2.original$beta[,scad.full.removal.part2.original$min][-1])

var.scad.removal.part2.original <- which(coeff.scad.removal.part2.original!=0)

################################################################################ 
# MCP 

#full
mcp.full.original <- cv.ncvreg(as.matrix(X), y, family="gaussian", penalty="MCP") 

mcp.full.obj.original <- mcp.full.original$fit 

coeff.mcp.full.original <- as.numeric(mcp.full.obj.original$beta[,mcp.full.original$min][-1])

var.mcp.full.original <- which(coeff.mcp.full.original!=0)

#removal part1 
y.mcp.removal.part1.original <- y[-ClusMIP_MCP]
X.mcp.removal.part1.original <- X[-ClusMIP_MCP,]

mcp.full.removal.part1.original <- cv.ncvreg(as.matrix(X.mcp.removal.part1.original), y.mcp.removal.part1.original, family="gaussian", penalty="MCP") 

mcp.full.obj.removal.part1.original <- mcp.full.removal.part1.original$fit 

coeff.mcp.removal.part1.original <- as.numeric(mcp.full.obj.removal.part1.original$beta[,mcp.full.removal.part1.original$min][-1])

var.mcp.removal.part1.original <- which(coeff.mcp.removal.part1.original!=0)

#removal part2 
y.mcp.removal.part2.original <- y[-c(49,52,55)]
X.mcp.removal.part2.original <- X[-c(49,52,55),]

mcp.full.removal.part2.original <- cv.ncvreg(as.matrix(X.mcp.removal.part2.original), y.mcp.removal.part2.original, family="gaussian", penalty="MCP") 

mcp.full.obj.removal.part2.original <- mcp.full.removal.part2.original$fit 

coeff.mcp.removal.part2.original <- as.numeric(mcp.full.obj.removal.part2.original$beta[,mcp.full.removal.part2.original$min][-1])

var.mcp.removal.part2.original <- which(coeff.mcp.removal.part2.original!=0)

################################################################################ 
#DF(LASSO)

#removal part 1 

y.df.removal.part1 <- as.vector(dat$Y)[-DF]

X.df.removal.part1 <- dat[-DF,-1]

df.removal.part1 <- cv.glmnet(as.matrix(X.df.removal.part1), y.df.removal.part1, alpha=1, family="gaussian", nfolds=10)   

coeff.df.removal.part1 <- as.numeric(coef(df.removal.part1, s="lambda.min")[-1]) 

var.df.removal.part1 <- which(coeff.df.removal.part1!=0)

y.hat.df.removal.part1 <- as.matrix(X.df.removal.part1) %*% coeff.df.removal.part1

pred.error.df.removal.part1 <- sum((y.df.removal.part1-y.hat.df.removal.part1)^2)

pred.error.df.removal.part1

cor(y.df.removal.part1, y.hat.df.removal.part1)




################################################################################ 
# ENET

# full 
coeff.enet.full.original <- enet.fit(X=as.matrix(X), Y=y, n_folds=10, dat.family="gaussian")$enet.coef   

var.enet.full.original <- which(coeff.enet.full.original!=0)

# removal part1 
y.enet.removal.part1.original <- y[-ClusMIP_enet]
X.enet.removal.part1.original <- X[-ClusMIP_enet,]

coeff.enet.removal.part1.original <- enet.fit(X=as.matrix(X.enet.removal.part1.original), Y=y.enet.removal.part1.original, n_folds=10, dat.family="gaussian")$enet.coef   

var.enet.removal.part1.original <- which(coeff.enet.removal.part1.original!=0)

#removal part2 
y.enet.removal.part2.original <- y[-c(49,52,55)]
X.enet.removal.part2.original <- X[-c(49,52,55),]

coeff.enet.removal.part2.original <- enet.fit(X=as.matrix(X.enet.removal.part2.original), Y=y.enet.removal.part2.original, n_folds=10, dat.family="gaussian")$enet.coef   

var.enet.removal.part2.original <- which(coeff.enet.removal.part2.original!=0)

################################################################################ 

# Effect on model selection

# ClusMIP(LASSO)
Name <- rep(c("Full", "Remove_Individual", "Remove_Common"),each=ncol(X))

var <- rep(1:ncol(X),3)

full <- rep(0,ncol(X))

b1 <- full
b1[var.lasso.full.original] <- 1

b2 <- full
b2[var.lasso.removal.part1.original] <- 1

b3 <- full
b3[var.lasso.removal.part2.original] <- 1

selection <- c(b1,b2,b3)

mat.lasso.selection <- data.frame(Full=b1, ReducedI=b2, ReducedII=b3)
write.csv(mat.lasso.selection,"selection_lasso.csv")

mat.lasso.coef <- data.frame(Full=coeff.lasso.full.original, 
                             ReducedI=coeff.lasso.removal.part1.original, 
                             ReducedII=coeff.lasso.removal.part2.original)
write.csv(mat.lasso.coef, "coef_lasso.csv")

dat.ClusMIP.lasso.selection <- data.frame(Procedure=Name, Index=var, Selection=selection)
dat.ClusMIP.lasso.selection$Procedure <- factor(dat.ClusMIP.lasso.selection$Procedure, levels = c("Remove_Common", "Remove_Individual", "Full"))

#ClusMIP(ENET)
Name <- rep(c("Full", "Remove_Individual", "Remove_Common"),each=ncol(X))

var <- rep(1:ncol(X),3)

full <- rep(0,ncol(X))

b1 <- full
b1[var.enet.full.original] <- 1

b2 <- full
b2[var.enet.removal.part1.original] <- 1

b3 <- full
b3[var.enet.removal.part2.original] <- 1

selection <- c(b1,b2,b3)

mat.enet.selection <- data.frame(Full=b1, ReducedI=b2, ReducedII=b3)
write.csv(mat.enet.selection,"selection_enet.csv")

mat.enet.coef <- data.frame(Full=coeff.enet.full.original, 
                             ReducedI=coeff.enet.removal.part1.original, 
                             ReducedII=coeff.enet.removal.part2.original)
write.csv(mat.enet.coef, "coef_enet.csv")

dat.ClusMIP.enet.selection <- data.frame(Procedure=Name, Index=var, Selection=selection)
dat.ClusMIP.enet.selection$Procedure <- factor(dat.ClusMIP.enet.selection$Procedure, levels = c("Remove_Common", "Remove_Individual", "Full"))

#ClusMIP(SCAD)
Name <- rep(c("Full", "Remove_Individual", "Remove_Common"),each=ncol(X))

var <- rep(1:ncol(X),3)

full <- rep(0,ncol(X))

b1 <- full
b1[var.scad.full.original] <- 1

b2 <- full
b2[var.scad.removal.part1.original] <- 1

b3 <- full
b3[var.scad.removal.part2.original] <- 1

selection <- c(b1,b2,b3)

mat.scad.selection <- data.frame(Full=b1, ReducedI=b2, ReducedII=b3)
write.csv(mat.scad.selection,"selection_scad.csv")

mat.scad.coef <- data.frame(Full=coeff.scad.full.original, 
                            ReducedI=coeff.scad.removal.part1.original, 
                            ReducedII=coeff.scad.removal.part2.original)
write.csv(mat.scad.coef, "coef_scad.csv")

dat.ClusMIP.scad.selection <- data.frame(Procedure=Name, Index=var, Selection=selection)
dat.ClusMIP.scad.selection$Procedure <- factor(dat.ClusMIP.scad.selection$Procedure, levels = c("Remove_Common", "Remove_Individual", "Full"))

#ClusMIP(MCP)
Name <- rep(c("Full", "Remove_Individual", "Remove_Common"),each=ncol(X))

var <- rep(1:ncol(X),3)

full <- rep(0,ncol(X))

b1 <- full
b1[var.mcp.full.original] <- 1

b2 <- full
b2[var.mcp.removal.part1.original] <- 1

b3 <- full
b3[var.mcp.removal.part2.original] <- 1

selection <- c(b1,b2,b3)

mat.mcp.selection <- data.frame(Full=b1, ReducedI=b2, ReducedII=b3)
write.csv(mat.mcp.selection,"selection_mcp.csv")

mat.mcp.coef <- data.frame(Full=coeff.mcp.full.original, 
                            ReducedI=coeff.mcp.removal.part1.original, 
                            ReducedII=coeff.mcp.removal.part2.original)
write.csv(mat.mcp.coef, "coef_mcp.csv")

dat.ClusMIP.mcp.selection <- data.frame(Procedure=Name, Index=var, Selection=selection)
dat.ClusMIP.mcp.selection$Procedure <- factor(dat.ClusMIP.mcp.selection$Procedure, levels = c("Remove_Common", "Remove_Individual", "Full"))

#Plot 
jpeg(filename=paste0("RealDataI_Selection.jpeg"), width=2400, height=1600) 

fig1 <- ggplot(dat.ClusMIP.lasso.selection, aes(x=Index, y=Procedure, fill=Selection)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") + 
  theme(axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.text.x = element_text(color="black",  size=20),
        axis.text.y = element_text(color="black",  size=20),
        legend.title = element_text(size=20),
        #legend.title = element_blank(),  
        legend.text = element_text(size=20),
        panel.border = element_rect(color = "black", fill = NA, size = 1))#,
        #plot.margin = unit(c(3,3,3,3), "lines"))  

fig2 <- ggplot(dat.ClusMIP.enet.selection, aes(x=Index, y=Procedure, fill=Selection)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") + 
  theme(axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.text.x = element_text(color="black",  size=20),
        axis.text.y = element_text(color="black",  size=20),
        legend.title = element_text(size=20),
        #legend.title = element_blank(),  
        legend.text = element_text(size=20),
        panel.border = element_rect(color = "black", fill = NA, size = 1))#,
#plot.margin = unit(c(3,3,3,3), "lines")) 

fig3 <- ggplot(dat.ClusMIP.scad.selection, aes(x=Index, y=Procedure, fill=Selection)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") + 
  theme(axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.text.x = element_text(color="black",  size=20),
        axis.text.y = element_text(color="black",  size=20),
        legend.title = element_text(size=20),
        #legend.title = element_blank(),  
        legend.text = element_text(size=20),
        panel.border = element_rect(color = "black", fill = NA, size = 1))#,
#plot.margin = unit(c(3,3,3,3), "lines")) 


fig4 <- ggplot(dat.ClusMIP.mcp.selection, aes(x=Index, y=Procedure, fill=Selection)) + 
  geom_tile() +
  scale_fill_gradient(low="white", high="red") + 
  theme(axis.title.x = element_text(color="black", size=20, face="bold"),
        axis.title.y = element_text(color="black", size=20, face="bold"),
        axis.text.x = element_text(color="black",  size=20),
        axis.text.y = element_text(color="black",  size=20),
        legend.title = element_text(size=20),
        #legend.title = element_blank(),  
        legend.text = element_text(size=20),
        panel.border = element_rect(color = "black", fill = NA, size = 1))#,
#plot.margin = unit(c(3,3,3,3), "lines")) 

final <- ggarrange(fig1, fig2, fig3, fig4,   
                   labels = c("ClusMIP(LASSO)", "ClusMIP(ENET)", "ClusMIP(SCAD)", "ClusMIP(MCP)"),
                   ncol = 2, nrow = 2, vjust=-2, common.legend = TRUE, font.label=list(color="black",size=22))

#legend="bottom",

annotate_figure(final, top = text_grob("Model Selection Outcomes",
                                       color = "blue", face = "bold", size = 26))

dev.off()

################################################################################

#cross validation and prediction correlation 

switch <- c(1,0,1,1)

cv.fun <- function(folds,B,X,y,ClusMIP_LASSO,ClusMIP_enet, ClusMIP_SCAD, ClusMIP_MCP, common, switch, output.dir){

  n <- nrow(X)
  
  cor.mat <- matrix(0,nrow=4,ncol=3)
  
  for(i in 1:B){
    
    test.index <- sort(sample(1:n, size=floor(n/5), replace = FALSE, prob = NULL), decreasing=FALSE)
    
    train.index <- setdiff(1:n,test.index)
    
    X.full.train <- X[train.index,]
    y.full.train <- y[train.index]
    
    X.full.test <- X[test.index,]
    y.full.test <- y[test.index]
    
    #write.csv(test.index, paste0(output.dir,"test_index.csv"))
    
    #write.csv(train.index, paste0(output.dir,"train_index.csv"))
    
    #write.csv(X.full.train, paste0(output.dir,"X_full_train.csv"))
    
    #write.csv(y.full.train, paste0(output.dir,"y_full_train.csv"))
    
    #write.csv(X.full.test, paste0(output.dir,"X_full_test.csv"))
    
    #write.csv(y.full.test, paste0(output.dir,"y_full_test.csv"))
    
    
    ############################################################################
    if(switch[1]==1){
    # LASSO: full 
    lasso.full <- cv.glmnet(as.matrix(X.full.train), y.full.train, alpha=1, family="gaussian", nfolds=10)   
    coeff.lasso.full <- as.numeric(coef(lasso.full, s="lambda.min")[-1]) 
    #var.lasso.full <- which(coeff.lasso.full!=0)
    y.hat.lasso.full <- as.matrix(X.full.test) %*% coeff.lasso.full
    cor.mat[1,1] <- cor.mat[1,1] + as.numeric(cor(y.full.test, y.hat.lasso.full))
    
    #write.csv(y.hat.lasso.full, paste0(output.dir,"y_hat_lasso_full.csv"))
    
    #write.csv(coeff.lasso.full, paste0(output.dir,"coeff_lasso_full.csv"))
    
    # LASSO: removal part1 
    index.lasso.removal.part1 <- setdiff(1:n,ClusMIP_LASSO)
    
    lasso.removal.part1.test.index <- sort(sample(index.lasso.removal.part1, 
                                                  size=floor(length(index.lasso.removal.part1)/5), replace = FALSE), 
                                                  decreasing=FALSE)
    
    lasso.removal.part1.train.index <- setdiff(index.lasso.removal.part1,lasso.removal.part1.test.index)
    
    #write.csv(lasso.removal.part1.test.index, paste0(output.dir,"lasso_removal_part1_test_index.csv"))
    
    X.lasso.removal.part1.train <- X[lasso.removal.part1.train.index,]
    y.lasso.removal.part1.train <- y[lasso.removal.part1.train.index]
    
    X.lasso.removal.part1.test <- X[lasso.removal.part1.test.index,]
    y.lasso.removal.part1.test <- y[lasso.removal.part1.test.index]
    
    lasso.removal.part1 <- cv.glmnet(as.matrix(X.lasso.removal.part1.train), y.lasso.removal.part1.train, alpha=1, family="gaussian", nfolds=10)   
    coeff.lasso.removal.part1 <- as.numeric(coef(lasso.removal.part1, s="lambda.min")[-1]) 
    #var.lasso.removal.part1 <- which(coeff.lasso.removal.part1!=0)
    y.hat.lasso.removal.part1 <- as.matrix(X.lasso.removal.part1.test) %*% coeff.lasso.removal.part1
    cor.mat[1,2] <- cor.mat[1,2] + cor(y.lasso.removal.part1.test, y.hat.lasso.removal.part1)
    
    #write.csv(y.hat.lasso.removal.part1, paste0(output.dir,"y_hat_lasso_removal_part1.csv"))
    
    #write.csv(coeff.lasso.removal.part1, paste0(output.dir,"coeff_lasso_removal_part1.csv"))
    
    # LASSO: removal part2 
    index.lasso.removal.part2 <- setdiff(1:n,common)
    
    lasso.removal.part2.test.index <- sort(sample(index.lasso.removal.part2, 
                                                  size=floor(length(index.lasso.removal.part2)/5), replace = FALSE), 
                                                  decreasing=FALSE)
    
    lasso.removal.part2.train.index <- setdiff(index.lasso.removal.part2,lasso.removal.part2.test.index)
    
    #write.csv(lasso.removal.part2.test.index, paste0(output.dir,"lasso_removal_part2_test_index.csv"))
    
    X.lasso.removal.part2.train <- X[lasso.removal.part2.train.index,]
    y.lasso.removal.part2.train <- y[lasso.removal.part2.train.index]
    
    X.lasso.removal.part2.test <- X[lasso.removal.part2.test.index,]
    y.lasso.removal.part2.test <- y[lasso.removal.part2.test.index]
    
    lasso.removal.part2 <- cv.glmnet(as.matrix(X.lasso.removal.part2.train), y.lasso.removal.part2.train, alpha=1, family="gaussian", nfolds=10)   
    coeff.lasso.removal.part2 <- as.numeric(coef(lasso.removal.part2, s="lambda.min")[-1]) 
    var.lasso.removal.part2 <- which(coeff.lasso.removal.part2!=0)
    y.hat.lasso.removal.part2 <- as.matrix(X.lasso.removal.part2.test) %*% coeff.lasso.removal.part2
    cor.mat[1,3] <- cor.mat[1,3] + cor(y.lasso.removal.part2.test, y.hat.lasso.removal.part2)
    
    #write.csv(y.hat.lasso.removal.part2, paste0(output.dir,"y_hat_lasso_removal_part2.csv"))
    
    #write.csv(coeff.lasso.removal.part2, paste0(output.dir,"coeff_lasso_removal_part2.csv"))
    
    }
    
    ############################################################################
    # ENET: full 
    if(switch[2]==1){
    coeff.enet.full <- enet.fit(X=as.matrix(X.full.train), Y=y.full.train, n_folds=10, dat.family="gaussian")$enet.coef 
    var.enet.full <- which(coeff.enet.full!=0)
    y.hat.enet.full <- as.matrix(X.full.test) %*% coeff.enet.full
    cor.mat[2,1] <- cor.mat[2,1] + cor(y.full.test, y.hat.enet.full)
    
    #write.csv(y.hat.enet.full, paste0(output.dir,"y_hat_enet_full.csv"))
    
    #write.csv(coeff.enet.full, paste0(output.dir,"coeff_enet_full.csv"))
    
    # ENET: removal part1 
    index.enet.removal.part1 <- setdiff(1:n,ClusMIP_enet)
    
    enet.removal.part1.test.index <- sort(sample(index.enet.removal.part1, 
                                          size=floor(length(index.enet.removal.part1)/5), replace = FALSE), 
                                          decreasing=FALSE)
    
    enet.removal.part1.train.index <- setdiff(index.enet.removal.part1,enet.removal.part1.test.index)
    
    #write.csv(enet.removal.part1.test.index, paste0(output.dir,"enet_removal_part1_test_index.csv"))
    
    X.enet.removal.part1.train <- X[enet.removal.part1.train.index,]
    y.enet.removal.part1.train <- y[enet.removal.part1.train.index]
    
    X.enet.removal.part1.test <- X[enet.removal.part1.test.index,]
    y.enet.removal.part1.test <- y[enet.removal.part1.test.index]
    
    coeff.enet.removal.part1 <- enet.fit(X=as.matrix(X.enet.removal.part1.train), Y=y.enet.removal.part1.train, n_folds=10, dat.family="gaussian")$enet.coef   
    var.enet.removal.part1 <- which(coeff.enet.removal.part1!=0)
    y.hat.enet.removal.part1 <- as.matrix(X.enet.removal.part1.test) %*% coeff.enet.removal.part1
    cor.mat[2,2] <- cor.mat[2,2] + cor(y.enet.removal.part1.test, y.hat.enet.removal.part1)
    
    #write.csv(y.hat.enet.removal.part1, paste0(output.dir,"y_hat_enet_removal_part1.csv"))
    
    #write.csv(coeff.enet.removal.part1, paste0(output.dir,"coeff_enet_removal_part1.csv"))
    
    #ENET: removal part2 
    index.enet.removal.part2 <- setdiff(1:n,common)
    
    enet.removal.part2.test.index <- sort(sample(index.enet.removal.part2, 
                                          size=floor(length(index.enet.removal.part2)/5), replace = FALSE), 
                                          decreasing=FALSE)
    
    enet.removal.part2.train.index <- setdiff(index.enet.removal.part2,enet.removal.part2.test.index)
    
    #write.csv(enet.removal.part2.test.index, paste0(output.dir,"enet_removal_part2_test_index.csv"))
    
    X.enet.removal.part2.train <- X[enet.removal.part2.train.index,]
    y.enet.removal.part2.train <- y[enet.removal.part2.train.index]
    
    X.enet.removal.part2.test <- X[enet.removal.part2.test.index,]
    y.enet.removal.part2.test <- y[enet.removal.part2.test.index]
    
    coeff.enet.removal.part2 <- enet.fit(X=as.matrix(X.enet.removal.part2.train), Y=y.enet.removal.part2.train, n_folds=10, dat.family="gaussian")$enet.coef   
    var.enet.removal.part2 <- which(coeff.enet.removal.part2!=0)
    y.hat.enet.removal.part2 <- as.matrix(X.enet.removal.part2.test) %*% coeff.enet.removal.part2
    cor.mat[2,3] <- cor.mat[2,3] + cor(y.enet.removal.part2.test, y.hat.enet.removal.part2)
    
    #write.csv(y.hat.enet.removal.part2, paste0(output.dir,"y_hat_enet_removal_part2.csv"))
    
    #write.csv(coeff.enet.removal.part2, paste0(output.dir,"coeff_enet_removal_part2.csv"))
    }
    
    ############################################################################
    # SCAD: full
    if(switch[3]==1){
    scad.full <- cv.ncvreg(as.matrix(X.full.train), y.full.train, family="gaussian", penalty="SCAD") 
    scad.full.obj <- scad.full$fit 
    coeff.scad.full <- as.numeric(scad.full.obj$beta[,scad.full$min][-1])
    var.scad.full <- which(coeff.scad.full!=0)
    y.hat.scad.full <- as.matrix(X.full.test) %*% coeff.scad.full
    cor.mat[3,1] <- cor.mat[3,1] + cor(y.full.test,y.hat.scad.full)
    
    #write.csv(y.hat.scad.full, paste0(output.dir,"y_hat_scad_full.csv"))
    
    #write.csv(coeff.scad.full, paste0(output.dir,"coeff_scad_full.csv"))
    
    # SCAD: removal part1
    index.scad.removal.part1 <- setdiff(1:n,ClusMIP_SCAD)
    
    scad.removal.part1.test.index <- sort(sample(index.scad.removal.part1, 
                                                 size=floor(length(index.scad.removal.part1)/5), replace = FALSE), 
                                          decreasing=FALSE)
    
    scad.removal.part1.train.index <- setdiff(index.scad.removal.part1,scad.removal.part1.test.index)
    
    #write.csv(scad.removal.part1.test.index, paste0(output.dir,"scad_removal_part1_test_index.csv"))
    
    X.scad.removal.part1.train <- X[scad.removal.part1.train.index,]
    y.scad.removal.part1.train <- y[scad.removal.part1.train.index]
    
    X.scad.removal.part1.test <- X[scad.removal.part1.test.index,]
    y.scad.removal.part1.test <- y[scad.removal.part1.test.index]
    
    scad.full.removal.part1 <- cv.ncvreg(as.matrix(X.scad.removal.part1.train), 
                                         y.scad.removal.part1.train, family ="gaussian", penalty="SCAD") 
    scad.full.obj.removal.part1 <- scad.full.removal.part1$fit 
    coeff.scad.removal.part1 <- as.numeric(scad.full.obj.removal.part1$beta[,scad.full.removal.part1$min][-1])
    var.scad.removal.part1 <- which(coeff.scad.removal.part1!=0)
    y.hat.scad.removal.part1 <- as.matrix(X.scad.removal.part1.test) %*% coeff.scad.removal.part1
    cor.mat[3,2] <- cor.mat[3,2] +cor(y.scad.removal.part1.test, y.hat.scad.removal.part1)
    
    #write.csv(y.hat.scad.removal.part1, paste0(output.dir,"y_hat_scad_removal_part1.csv"))
    
    #write.csv(coeff.scad.removal.part1, paste0(output.dir,"coeff_scad_removal_part1.csv"))
    
    # SCAD: removal part2
    index.scad.removal.part2 <- setdiff(1:n,common)
    
    scad.removal.part2.test.index <- sort(sample(index.scad.removal.part2, 
                                                 size=floor(length(index.scad.removal.part2)/5), replace = FALSE), 
                                          decreasing=FALSE)
    
    scad.removal.part2.train.index <- setdiff(index.scad.removal.part2,scad.removal.part2.test.index)
    
    #write.csv(scad.removal.part2.test.index, paste0(output.dir,"scad_removal_part2_test_index.csv"))
    
    X.scad.removal.part2.train <- X[scad.removal.part2.train.index,]
    y.scad.removal.part2.train <- y[scad.removal.part2.train.index]
    
    X.scad.removal.part2.test <- X[scad.removal.part2.test.index,]
    y.scad.removal.part2.test <- y[scad.removal.part2.test.index]
    
    scad.full.removal.part2 <- cv.ncvreg(as.matrix(X.scad.removal.part2.train), 
                                         y.scad.removal.part2.train, family ="gaussian", penalty="SCAD") 
    scad.full.obj.removal.part2 <- scad.full.removal.part2$fit 
    coeff.scad.removal.part2 <- as.numeric(scad.full.obj.removal.part2$beta[,scad.full.removal.part2$min][-1])
    var.scad.removal.part2 <- which(coeff.scad.removal.part2!=0)
    y.hat.scad.removal.part2 <- as.matrix(X.scad.removal.part2.test) %*% coeff.scad.removal.part2
    cor.mat[3,3] <- cor.mat[3,3] +cor(y.scad.removal.part2.test, y.hat.scad.removal.part2)
    
    #write.csv(y.hat.scad.removal.part2, paste0(output.dir,"y_hat_scad_removal_part2.csv"))
    
    #write.csv(coeff.scad.removal.part2, paste0(output.dir,"coeff_scad_removal_part2.csv"))
    
    }
    
    ############################################################################
    # MCP: full
    if(switch[4]==1){
    mcp.full <- cv.ncvreg(as.matrix(X.full.train), y.full.train, family="gaussian", penalty="MCP") 
    mcp.full.obj <- mcp.full$fit 
    coeff.mcp.full <- as.numeric(mcp.full.obj$beta[,mcp.full$min][-1])
    var.mcp.full <- which(coeff.mcp.full!=0)
    y.hat.mcp.full <- as.matrix(X.full.test) %*% coeff.mcp.full
    cor.mat[4,1] <- cor.mat[4,1] + cor(y.full.test,y.hat.mcp.full)
    
    #write.csv(y.hat.mcp.full, paste0(output.dir,"y_hat_mcp_full.csv"))
    
    #write.csv(coeff.mcp.full, paste0(output.dir,"coeff_mcp_full.csv"))
    
    #MCP: removal part1 
    index.mcp.removal.part1 <- setdiff(1:n,ClusMIP_MCP)
    
    mcp.removal.part1.test.index <- sort(sample(index.mcp.removal.part1, 
                                                 size=floor(length(index.mcp.removal.part1)/5), replace = FALSE), 
                                          decreasing=FALSE)
    
    mcp.removal.part1.train.index <- setdiff(index.mcp.removal.part1,mcp.removal.part1.test.index)
    
    #write.csv(mcp.removal.part1.test.index, paste0(output.dir,"mcp_removal_part1_test_index.csv"))
    
    X.mcp.removal.part1.train <- X[mcp.removal.part1.train.index,]
    y.mcp.removal.part1.train <- y[mcp.removal.part1.train.index]
    
    X.mcp.removal.part1.test <- X[mcp.removal.part1.test.index,]
    y.mcp.removal.part1.test <- y[mcp.removal.part1.test.index]
    
    mcp.full.removal.part1 <- cv.ncvreg(as.matrix(X.mcp.removal.part1.train), 
                                        y.mcp.removal.part1.train, family="gaussian", penalty="MCP") 
    mcp.full.obj.removal.part1 <- mcp.full.removal.part1$fit 
    coeff.mcp.removal.part1 <- as.numeric(mcp.full.obj.removal.part1$beta[,mcp.full.removal.part1$min][-1])
    var.mcp.removal.part1 <- which(coeff.mcp.removal.part1!=0)
    y.hat.mcp.removal.part1 <- as.matrix(X.mcp.removal.part1.test) %*% coeff.mcp.removal.part1
    cor.mat[4,2] <- cor.mat[4,2] + cor(y.mcp.removal.part1.test, y.hat.mcp.removal.part1)
    
    #write.csv(y.hat.mcp.removal.part1, paste0(output.dir,"y_hat_mcp_removal_part1.csv"))
    
    #write.csv(coeff.mcp.removal.part1, paste0(output.dir,"coeff_mcp_removal_part1.csv"))
    
    #MCP: removal part2 
    index.mcp.removal.part2 <- setdiff(1:n,common)
    
    mcp.removal.part2.test.index <- sort(sample(index.mcp.removal.part2, 
                                                size=floor(length(index.mcp.removal.part2)/5), replace = FALSE), 
                                         decreasing=FALSE)
    
    mcp.removal.part2.train.index <- setdiff(index.mcp.removal.part2,mcp.removal.part2.test.index).
    
    #write.csv(mcp.removal.part2.test.index, paste0(output.dir,"mcp_removal_part2_test_index.csv"))
    
    X.mcp.removal.part2.train <- X[mcp.removal.part2.train.index,]
    y.mcp.removal.part2.train <- y[mcp.removal.part2.train.index]
    
    X.mcp.removal.part2.test <- X[mcp.removal.part2.test.index,]
    y.mcp.removal.part2.test <- y[mcp.removal.part2.test.index]
    
    mcp.full.removal.part2 <- cv.ncvreg(as.matrix(X.mcp.removal.part2.train), 
                                        y.mcp.removal.part2.train, family="gaussian", penalty="MCP") 
    mcp.full.obj.removal.part2 <- mcp.full.removal.part2$fit 
    coeff.mcp.removal.part2 <- as.numeric(mcp.full.obj.removal.part2$beta[,mcp.full.removal.part2$min][-1])
    var.mcp.removal.part2 <- which(coeff.mcp.removal.part2!=0)
    y.hat.mcp.removal.part2 <- as.matrix(X.mcp.removal.part2.test) %*% coeff.mcp.removal.part2
    cor.mat[4,3] <- cor.mat[4,3] + cor(y.mcp.removal.part2.test, y.hat.mcp.removal.part2)
    
    #write.csv(y.hat.mcp.removal.part2, paste0(output.dir,"y_hat_mcp_removal_part2.csv"))
    
    #write.csv(coeff.mcp.removal.part2, paste0(output.dir,"coeff_mcp_removal_part2.csv"))
    
    }
    
  }
  
  cor.mat <- cor.mat/B
  
  rownames(cor.mat) <- c("ClusMIP(LASSO)", "ClusMIP(ENET)", "ClusMIP(SCAD)", "ClusMIP(MCP)")
  colnames(cor.mat) <- c("Full", "Removal1", "Removal2")
  
  write.csv(cor.mat, paste0(output.dir,"cor_mat.csv"))
  
  returnList <- list("cor.mat"=cor.mat)
  
  return(returnList)
  
}

Method <- rep(c("ClusMIP(LASSO)", "ClusMIP(ENET)", "ClusMIP(SCAD)", "ClusMIP(MCP)"), each=3)
Operation <- rep(c("Full", "Removal1", "Removal2"), 4)
Correlation <- as.numeric(t(cor.mat))
df.prediction <- data.frame(Method=Method, Operation=Operation, Correlation=Correlation)
df.prediction$Method <- factor(df.prediction$Method, levels = c("ClusMIP(LASSO)", "ClusMIP(ENET)", "ClusMIP(SCAD)", "ClusMIP(MCP)"))
df.prediction$Operation <- factor(df.prediction$Operation, levels = c("Full", "Removal1", "Removal2"))

#Plot 
jpeg(filename=paste0(output.dir,"RealDataI_Prediction_Correlation.jpeg"), width=2400, height=1600) 

txt.ylab <- expression(paste("Correlation between Y and ", hat(Y))) 
#txt.xlab.removal1 <- expression(paste("Removal of ", hat(I)[ClusMIP(Selector)]^'infl')) 
#txt.xlab.removal2 <- expression(paste("Removal of {49,52,55}")) 

txt.xlab.removal1 <- expression(paste("Reduced Dataset\nupon removing {", Z[r], ": r in ", hat(I)[ClusMIP(Selector)]^'infl', "}")) 
txt.xlab.removal2 <- expression(paste("Removal of {49,52,55}")) 

fig1 <- ggplot(df.prediction, aes(x=Operation, y=Correlation, group=Method)) + 
        geom_line(linetype = "dashed") + 
        geom_point(aes(color=Method), size = 6)  +
        xlab("") + ylab(txt.ylab) + 
        #geom_tile() +
        #scale_fill_gradient(low="white", high="red") + 
        theme(axis.title.x = element_text(color="black", size=36, face="bold"),
              axis.title.y = element_text(color="black", size=36, face="bold"),
              axis.text.x = element_text(color="black",  size=36, vjust=-0.5),
              axis.text.y = element_text(color="black",  size=36),
              legend.title = element_text(size=36,face="bold"),
        #legend.title = element_blank(),  
        legend.text = element_text(size=36),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = unit(c(3,3,3,3), "lines"))  + 
        scale_y_continuous(breaks = seq(0.5, 1.0, 0.05)) + 
        scale_colour_discrete(name="Detection\nProcedures", labels=c("ClusMIP(LASSO)", "ClusMIP(ENET)", "ClusMIP(SCAD)", "ClusMIP(MCP)")) +
        scale_x_discrete(labels=c("Full Dataset", txt.xlab.removal1, txt.xlab.removal2))

print(fig1)

dev.off()

Method <- rep(c("ClusMIP(LASSO)", "ClusMIP(ENET)", "ClusMIP(SCAD)", "ClusMIP(MCP)"), each=3)
Operation <- rep(c("Full", "Removal1", "Removal2"), 4)
Correlation <- as.numeric(t(cor.mat))
df.prediction <- data.frame(Method=Method, Operation=Operation, Correlation=Correlation)
df.prediction$Method <- factor(df.prediction$Method, levels = c("ClusMIP(LASSO)", "ClusMIP(ENET)", "ClusMIP(SCAD)", "ClusMIP(MCP)"))
df.prediction$Operation <- factor(df.prediction$Operation, levels = c("Full", "Removal1", "Removal2"))

#Plot 
jpeg(filename=paste0(output.dir,"RealDataI_Prediction_Correlation.jpeg"), width=2400, height=1600) 

txt.ylab <- expression(paste("Correlation between Y and ", hat(Y))) 
#txt.xlab.removal1 <- expression(paste("Removal of ", hat(I)[ClusMIP(Selector)]^'infl')) 
#txt.xlab.removal2 <- expression(paste("Removal of {49,52,55}")) 

txt.xlab.removal1 <- expression(paste("Reduced Dataset\nupon removing {", Z[r], ": r in ", hat(I)[ClusMIP(Selector)]^'infl', "}")) 
txt.xlab.removal2 <- expression(paste("Removal of {49,52,55}")) 

fig1 <- ggplot(df.prediction, aes(x=Operation, y=Correlation, group=Method)) + 
  geom_line(linetype = "dashed") + 
  geom_point(aes(color=Method), size = 6)  +
  xlab("") + ylab(txt.ylab) + 
  #geom_tile() +
  #scale_fill_gradient(low="white", high="red") + 
  theme(axis.title.x = element_text(color="black", size=36, face="bold"),
        axis.title.y = element_text(color="black", size=36, face="bold"),
        axis.text.x = element_text(color="black",  size=36, vjust=-0.5),
        axis.text.y = element_text(color="black",  size=36),
        legend.title = element_text(size=36,face="bold"),
        #legend.title = element_blank(),  
        legend.text = element_text(size=36),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = unit(c(3,3,3,3), "lines"))  + 
  scale_y_continuous(breaks = seq(0.5, 1.0, 0.05)) + 
  scale_colour_discrete(name="Detection\nProcedures", labels=c("ClusMIP(LASSO)", "ClusMIP(ENET)", "ClusMIP(SCAD)", "ClusMIP(MCP)")) +
  scale_x_discrete(labels=c("Full Dataset", txt.xlab.removal1, txt.xlab.removal2))

print(fig1)

dev.off()





