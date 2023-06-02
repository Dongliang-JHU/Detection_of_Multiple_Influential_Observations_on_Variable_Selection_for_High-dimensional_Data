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
# R Script Purpose: to compute the decision based on normal approximation 
#
# Reference: M.D.Judin, A Limit Theorem for Sums of Dependent Random Variables, Izv.Vyss Ucebn, Zaved Matematika, 
#            1962, no.1(26), 171-177, English Transl., Selected Transl.Math. Statist. and Probability,  Volume 5, 
#            American Mathematical Society, Providence, R.I., 1965, 308-314, MR 25 #1564
#
# Created on: October 04, 2021
#             
#
######################################################################################################################
# Function: normal.Decision 
#
# Purpose of Function: to compute the decision from normal approximation 
#
# Args:
#
#   dat             :  the count data vector (row vector)
#
#   scaling.factor  :  see Judin (1965), less than 1/4
#
#   alpha           :  type I error rate 
#
#
# Returns: 
#
#   decision.vec    :  a vector (length n) of 0 and 1 indicating flagging of influential points with 1 being flagged and 0 otherwise 
#

normal.Decision <- function(dat, scaling.factor, alpha){
  
  decision.vec <- rep(0,length(dat))
  
  dat.scale <- scaling.factor * dat
  
  upper <- mean(dat.scale) + qnorm(1-alpha/2) * sd(dat.scale)
  
  ind.infl <- which(dat.scale > upper)
  
  if( length(ind.infl) > 0){
    decision.vec[ind.infl] <- 1 
  }
    return(decision.vec)
}