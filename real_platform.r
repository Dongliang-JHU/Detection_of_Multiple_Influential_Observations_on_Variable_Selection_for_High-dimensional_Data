# Author: Dongliang Zhang (^)
#
# Supervisors: Professor Masoud Asgharian (*)
#              Professor Martin Lindquist (^)
#              Professor Mei-Cheng Wang (^)
#
# Current Affiliation: (*) Department of Mathematics and Statistics, McGill University, Montreal, Quebec, Canada 
#                      (^) Department of Biostatistics, Johns Hopkins University, Baltimore, Maryland, USA  
#
# Title of Project: Detection of Influential Observations on Variable Selection 
#
# R Script Purpose: Detection Algorithm for Real Data 
# 
# Created on  : August 24, 2022 
#
# Modified on : August 24, 2022 
#               August 25, 2022 
#
# Dependent R Scripts: dataGenerate.r, 
#                      enet_fit.r
#                      step1_multDetection.r,
#                      step1_multCluster.r,
#                      step2_filtering.r,
#                      normalDecision.r,
#                      cluster.r,
#                      binaryClassification.r,
#                      lasso_fit.r, 
#                      df_measure.r,
#                      testsEval.r; 
#
# R Packages :  MASS,
#               glmnet, 
#               scalreg,
#               ncvreg,
#               cluster, 
#               MIP (Zhao et al, 2019),
#               rlist (saving a list),
#               rje (is.subset),
#               tictoc (time), 
#               e1071(classification),
#               klaR (classification),
#               mpath 
#
##########################################################################################################################################
# Function: real.detection
#
# Args:
#   X                :  real design matrx 
#   y                :  real response vector  
#   detection.switch :  a switch control (0 or 1) for detection methods 
#                       Position 1  : MIP-VS + Clustering + LASSO 
#                       Position 2  : MIP-VS + Clustering + Scaled LASSO 
#                       Position 3  : MIP-VS + Clustering + Elastic Net  
#                       Position 4  : MIP-VS + Clustering + SCAD  
#                       Position 5  : MIP-VS + Clustering + MCP 
#                       Position 6  : MIP (Zhao et al, 2019)
#                       Position 7  : DF (Rajaratnam et al, 2019) 
#                       Position 8  : MIP-VS + Poisson + Clustering + LASSO + New Tuning Parameter (Jia and Xie, 2019)
#   alpha            :  FDR
#   n_folds          :  number of cross-validation 
#   n_subset         :  (MIP & MIP-VS) the number of subsets chosen at random to compute the Min and Max statistics, e.g n_subset=100
#   subset.size      :  (MIP & MIP-VS) the samples size in each subsets, e.g. subset_vol=n/2
#   K.vec            :  a vecotor of number of clusters to be determined, e.g. K.vec <- c(2,3)
#   family           :  gaussian or poisson 
#
# Returns: 
#   dec.mat          :  decision matrix (1: influential point, 0: non-influential point)
#   time.mat         :  time matrix 
#   infl.pos.cluster :  indices of suspected influential points 
#   infl.pos         :  list of indices of influential points detected 

real.detection <- function(X, y, cluster.type, detection.switch, alpha, K.vec, n_folds, n_subset, subset.size, family){
  
  n <- nrow(X)
  p <- ncol(X)
  
  Xy.dat <- as.data.frame(cbind(X, y))
  names(Xy.dat)<-c(paste("X", 1:p, sep=""),"Y")
  
  method.use.pos <- which(detection.switch==1)

  ################################################################################
  
  method.names <- c("LASSO", "SLASSO", "ENET", "SCAD", "MCP", "MIP", "DF", "LASSO_Poisson_New_Par")
  
  # times 
  time.mat <- matrix(NA, nrow=8, ncol=1)
  rownames(time.mat) <- method.names
  colnames(time.mat) <- "Time"
    
  # decision matrix 
  dec.mat <- matrix(0, nrow=8, ncol=n)
  rownames(dec.mat) <- method.names
  
  ################################################################################
  # clustering 
  
  if(cluster.type=="K.means"){
    clean.pos <- cluster.est(df=Xy.dat, K.vec) 
  } 
  
  ################################################################################
  # methods 
  
  # LASSO 
  if(detection.switch[1]==1){
    tic(paste0("MIP + VS + Clustering + LASSO"))
    obj.lasso <- mult.cluster(X=X, Y=y, n_folds, method.switch=c(1,0,0,0,0,0,0,0), dat.family=family, 
                              clean.pos, alpha)
    decision.lasso <- obj.lasso$dec.lasso
    print(which(decision.lasso==1))
    exectime <- toc()
    dec.mat[1,] <- decision.lasso
    time.mat[1,] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
    print(Sys.time())
  }
  
  # scaled LASSO
  if(detection.switch[2]==1){
    tic(paste0("MIP + VS + Clustering + Scaled LASSO"))
    obj.slasso <- mult.cluster(X=X, Y=y, n_folds, method.switch=c(0,0,1,0,0,0,0,0), dat.family=family, 
                               clean.pos, alpha)
    decision.slasso <- obj.slasso$dec.slasso
    print(which(decision.slasso==1))
    exectime <- toc()
    dec.mat[2,] <- decision.slasso 
    time.mat[2,] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
    print(Sys.time())
  }
    
  # elastic net 
  if(detection.switch[3]==1){
    tic(paste0("MIP + VS + Clustering + Elastic Net"))
    obj.enet <- mult.cluster(X=X, Y=y, n_folds, method.switch=c(0,0,0,1,0,0,0,0), dat.family=family, 
                             clean.pos, alpha)
    decision.enet <- obj.enet$dec.enet
    print(which(decision.enet==1))
    exectime <- toc()
    dec.mat[3,] <- decision.enet 
    time.mat[3,] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
    print(Sys.time())
  }
    
  # SCAD 
  if(detection.switch[4]==1){
    tic(paste0("MIP + VS + Clustering + SCAD"))
    obj.scad <- mult.cluster(X=X, Y=y, n_folds, method.switch=c(0,0,0,0,1,0,0,0), dat.family=family, 
                             clean.pos, alpha)
    decision.scad <- obj.scad$dec.scad
    print(which(decision.scad==1))
    exectime <- toc()
    dec.mat[4,] <- decision.scad 
    time.mat[4,] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
    print(Sys.time())
  }
    
  # MCP
  if(detection.switch[5]==1){
    tic(paste0("MIP + VS + Clustering + MCP"))
    obj.mcp <- mult.cluster(X=X, Y=y, n_folds, method.switch=c(0,0,0,0,0,1,0,0), dat.family=family, 
                            clean.pos, alpha)
    decision.mcp <- obj.mcp$dec.mcp
    print(which(decision.mcp==1))
    exectime <- toc()
    dec.mat[5,] <- decision.mcp
    time.mat[5,] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
    print(Sys.time())
  }
  
  # MIP (Zhao et al, 2019)
  if(detection.switch[6]==1){
    tic(paste0("MIP"))
    MIP.obj <- MIP(X=as.matrix(X), Y=t(as.matrix(y)), n, p, q=1, n_subset, subset.size, ep = 0.1, alpha)
    infl.pos.MIP <- sort(MIP.obj$inf_setfinal, decreasing=FALSE)
    print(infl.pos.MIP)
    exectime <- toc()
    dec.mat[6,infl.pos.MIP] <- 1  
    time.mat[6,] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
    print(Sys.time())
  }
  
  # DF (Rajaratnam et al, 2019) 
  if(detection.switch[7]==1){
    tic(paste0("DF"))
    decision.DF <- df.measure(X=X, Y=y)
    print(which(decision.DF==1))
    exectime <- toc()
    dec.mat[7,] <- decision.DF 
    time.mat[7,] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
    print(Sys.time())
  }
  
  # LASSO for Poisson with new tuning parameters (Jia and Xie, 2019)
  if(detection.switch[8]==1){
    tic(paste0("MIP + VS + Poisson + Clustering + LASSO + New Tuning Parameter"))
    obj.lasso.pois.new <- mult.cluster(X=X, Y=y, n_folds, method.switch=c(0,0,0,0,0,0,1,0), dat.family=family, 
                                       clean.pos, alpha)
    decision.lasso.pois.new <- obj.lasso.pois.new$dec.lasso.pois
    print(which(decision.lasso.pois.new==1))
    exectime <- toc()
    dec.mat[8,] <- decision.lasso.pois.new
    time.mat[8,] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
    print(Sys.time())
  }
  
  ################################################################################
  
  dec.mat.final <- dec.mat[method.use.pos,]
  time.mat.final <- time.mat[method.use.pos,]
  
  infl.pos.list <- lapply(1:nrow(dec.mat.final), function(i) which(dec.mat.final[i,]==1))
  
  returnList <- list("dec.mat"          = dec.mat.final,
                     "time.mat"         = time.mat.final, 
                     "infl.pos.cluster" = setdiff(1:n, clean.pos),
                     "infl.pos"         = infl.pos.list)
  
  return(returnList)  
  
}


