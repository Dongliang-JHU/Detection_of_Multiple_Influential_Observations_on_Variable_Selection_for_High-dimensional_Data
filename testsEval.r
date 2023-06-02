# Author: Dongliang Zhang (^)
#
# Supervisors: Professor Masoud Asgharian (*)
#              Professor Martin Lindquist (^)
#
# Current Affiliation: (1) Department of Mathematics and Statistics, McGill University, Montreal, Quebec, Canada (*)
#                      (2) Department of Biostatistics, Johns Hopkins University, Baltimore, Maryland, USA (^)
#
# Title of Project: Detection of Influential Cases on Model Selection 
#
# R Script Purpose: Power, False Positive Rate (FPR) and False Discovery Rate (FDR)
#
# Created on  : June 26, 2018
#
# Modified on : July 9, 2018
#               July 12, 2018
#               October 1, 2021
#               Ocotber 4, 2021
#
#########################################################################################################
#
# Function: testsEval()
#
# Purpose: to compute the power, false positive rate (FPR) and true positive rate (TPR); 
#
# Args:
#  mat.decision    :  a matrix of decisions ({0,1}'s) with dimension B times n; "1" indicates rejection;  
#  infl.index      :  a pre-supplied vector indicating indices of influential cases;
#  noninfl.index   :  a pre-supplied vector indicating indices of non-influential cases;
#
# Returns: 
#  power       :  power of a test;
#  FPR         :  false positive rate (FPR) of a test;
#  TPR         :  true positive rate (TPR) of a test; 
#
# Verification: 
#> mat.rejNonrej <- matrix(c(1,1,0,0,1,1,0,0,1,1,1,0),nrow=3)
#> mat.rejNonrej
#     [,1] [,2] [,3] [,4]
#[1,]    1    0    0    1
#[2,]    1    1    0    1
#[3,]    0    1    1    0
#> infl.index <- c(1,2)
#> B <- 3
#> n <- 4
#> power.row <- apply(mat.rejNonrej, 1, function(x) length(intersect(which(x!=0),infl.index)) ) 
#> power.row
#[1] 1 2 1
#> power.test <- sum(power.row) / (length(infl.index) * B);
#> power.test
#[1] 0.6666667


testsEval <- function(mat.decision, infl.index, noninfl.index){

 #H0: the case is not influential; 
 #H1: the case is influential; 
 #positive: rejecting H0;
 #negative: not rejecting H0; 
 
 #true positive(s); 
 TP <- apply(mat.decision, 1, function(x) length(intersect(which(x!=0), infl.index)) ); #x by default is each row of the matrix 
 
 #false positive(s); 
 FP <- apply(mat.decision, 1, function(x) length(intersect(which(x!=0), noninfl.index)) ); #x by default is each row of the matrix
 
 #true negative(s);
 TN <- apply(mat.decision, 1, function(x) length(intersect(which(x==0), noninfl.index)) ); #x by default is each row of the matrix
 
 #false negative(s);
 FN <- apply(mat.decision, 1, function(x) length(intersect(which(x==0), infl.index)) ); #x by default is each row of the matrix
 
 #total number of rejections 
 TR <- apply(mat.decision, 1, function(x) length(which(x!=0)) ); #x by default is each row of the matrix

 #power.test = total number of observations that are correctly rejected / total number of influential observations; 
 power.test <- sum(TP) / (length(infl.index) * nrow(mat.decision)); 
 
 #False Positive Rate (FPR) = FP / (FP + TN); 
 #probability of false alarms; 
 #1-Specificity
 FPR.test <- sum(FP) / (sum(FP) + sum(TN));
 
 #True Positive Rate (TPR) = TP / (TP + FN); 
 #probability of detection; 
 #a.k.a. Sensitivity 
 #same as the power of the test 
 #TPR.test <- sum(TP) / (sum(TP) + sum(FN));
 
 #False Discovery Rate 
 FDR <- sum(FP) / (sum(FP) + sum(TP))
   
 return(c(power.test, FPR.test, FDR)); 
 
}