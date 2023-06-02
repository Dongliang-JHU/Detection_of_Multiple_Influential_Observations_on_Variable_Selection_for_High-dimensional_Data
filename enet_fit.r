# Author: Dongliang Zhang (^)
#
# Supervisors: Professor Masoud Asgharian (*)
#              Professor Martin Lindquist (^)
#
# Current Affiliation: (1) Department of Mathematics and Statistics, McGill University, Montreal, Quebec, Canada (*)
#                      (2) Department of Biostatistics, Johns Hopkins University, Baltimore, Maryland, USA (^)
#
# Title of Project: Detection of Influential Cases on Variable Selection 
#
# R Script Purpose: fitting a elastic net penalized regression model to data 
#      
#
# Created on      : February 23, 2022
#
# Modified on     : February 24, 2022
#                   March 20, 2022
#                   March 28, 2022
#                   April 07, 2022 
#                   April 14, 2022 
#                   May 16, 2022 
#                   May 18, 2022 
#
# Dependent R Scripts: dataGenerate.r;  
#
# R Packages         : glmnet,
#                      ensr; 
#
#############################################################################################################################
# Function: enet.fit
#
# Args:
#
#   X                :  design matrx 
#   Y                :  response vector 
#   n_folds          :  number of cross-validation 
#
# Returns:  
#
#   

enet.fit <- function(X, Y, n_folds, dat.family){
  
  lambda.seq <- seq(0, 1, 0.1)
  
  search <- foreach(i=lambda.seq, .combine=rbind) %dopar% {
            cv.obj <- cv.glmnet(X, Y, family=dat.family, nfold=n_folds, type.measure="deviance", paralle=TRUE, alpha=i)
            data.frame(cvm=cv.obj$cvm[cv.obj$lambda==cv.obj$lambda.min], lambda.min=cv.obj$lambda.min, alpha=i)
            }
  
  tmp <- search[search$cvm==min(search$cvm), ]
  
  enet.obj <- glmnet(X, Y, family=dat.family, lambda=tmp$lambda.min, alpha=tmp$alpha)
  
  coeff <- as.numeric(coef(enet.obj))[-1] 
  
  returnList <- list("enet.coef"=coeff)
  
  return(returnList) 

}

#> search
#cvm lambda.min alpha
#1  52.68273  39.093353   0.0
#2  51.24968  10.144820   0.1
#3  50.51253   8.864180   0.2
#4  49.34370   6.190837   0.3
#5  48.35607   4.643127   0.4
#6  48.01050   3.714502   0.5
#7  48.41256   2.820430   0.6
#8  46.81909   2.779551   0.7
#9  46.07077   2.321564   0.8
#10 46.80542   1.713247   0.9
#11 46.12652   1.692258   1.0


