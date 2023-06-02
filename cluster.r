# Author: Dongliang Zhang (^)
#
# Supervisors: Professor Masoud Asgharian (*)
#              Professor Martin Lindquist (^)
#              Professor Mei-Cheng Wang (^)
#
# Current Affiliation: (1) Department of Mathematics and Statistics, McGill University, Montreal, Quebec, Canada (*)
#                      (2) Department of Biostatistics, Johns Hopkins University, Baltimore, Maryland, USA (^)
#
# Title of Project: Detection of Influential Observations on Variable Selection 
#
# R Script Purpose: Clustering   
# 
# Created on  : April 10, 2022 
#
# Modified on : April 10, 2022 
#               April 11, 2022
#               April 19, 2022 
#
# Dependent R Scripts: NA 
#
# R Packagess : stats; 
#
##########################################################################################################################################
# Function: cluster.est
#
# Reference:  
#
# Args:
#   df      :  data frame
#   K.vec   :  a vector of number of clusters to be performed, e.g. K.vec = c(2,3)
#    
# Returns: 
#   

cluster.est <- function(df, K.vec){
  
  #k: scalar 
  #k-means clustering
  #which cluster has the largest number of observations 
  myfunction <- function(df, k){
    
    clus.obj <- kmeans(df, k)
    
    cluster.pos <- clus.obj$cluster
    
    cluster.level <- unique(cluster.pos) 
    
    l.vec <- vector(length=length(cluster.level))
    
    for(i in 1:length(cluster.level)){
      l.vec[i] <- length(which(cluster.pos==cluster.level[i]))
    }
    
    pos <- which( cluster.pos== cluster.level[which.max(l.vec)] )
    
    clus.pos <- sort(pos, decreasing=FALSE)
    
    return(clus.pos)
  }
  
  pos.list <- lapply( K.vec, myfunction, df=df )
  clean.pos <- sort( Reduce(intersect, pos.list), decreasing=FALSE )

  return(clean.pos)
}

#Example: 
# a <- c(1,3,5,7,9)
# b <- c(3,6,8,9,10)
# c <- c(2,3,4,5,7,9)
# A <- list(a,b,c)
# Reduce(intersect,A)
#[1] 3 9

#> a <- c(4,5,2,2,2,2,2)
#> unique(a)
#[1] 4 5 2
