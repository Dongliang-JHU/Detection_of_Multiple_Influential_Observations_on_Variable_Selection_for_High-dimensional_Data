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
# R Script Purpose: Simulation Platform for Parallel Computing 
# 
# Created on  : May 04, 2023 
#
# Modified on : May 04, 2023 
#               May 05, 2023 
#               May 06, 2023
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
################################################################################

# Function: sim.platform.parallel
#
# Reference:  
#
# Args:
#   (DATA INPUT)
#   X                :  design matrx 
#   y                :  response vector 
#   family           :  "gaussian" or "poisson" 
#
#   (PERTURBATION)
#   perturb.type     :  perturbation type
#   y.infl.pos       :  indices of influential response 
#   x.infl.pos       :  indices of influential explanatory variables 
#   y.perturb        :  magnitude of perturbation on response vector  
#   x.perturb        :  magnitude of perturbation on explanatory variables  
#
#   (CONTROL SWITCH)
#   detection.switch :  a switch control (0 or 1) for detection methods 
#                       (linear models)
#                       Position 1  : MIP-VS + LASSO
#                       Position 2  : MIP-VS + Adaptive LASSO 
#                       Position 3  : MIP-VS + Scaled LASSO 
#                       Position 4  : MIP-VS + Elastic Net 
#                       Position 5  : MIP-VS + SCAD
#                       Position 6  : MIP-VS + MCP 
#                       (Clustering: linear models)
#                       Position 7  : MIP-VS + Clustering + LASSO 
#                       Position 8  : MIP-VS + Clustering + Adaptive LASSO 
#                       Position 9  : MIP-VS + Clustering + Scaled LASSO 
#                       Position 10 : MIP-VS + Clustering + Elastic Net  
#                       Position 11 : MIP-VS + Clustering + SCAD  
#                       Position 12 : MIP-VS + Clustering + MCP 
#                       Position 13 : MIP (Zhao et al, 2019)
#                       Position 14 : DF (Rajaratnam et al, 2019) 
#                       (Binary Classification)
#                       Position 15 : MIP-VS + Binary + LASSO 
#                       Position 16 : MIP-VS + Binary + SCAD
#                       Position 17 : MIP-VS + Binary + MCP
#                       (Binary Classification + Clustering)
#                       Position 18 : MIP-VS + Binary + Clustering + LASSO 
#                       Position 19 : MIP-VS + Binary + Clustering + SCAD
#                       Position 20 : MIP-VS + Binary + Clustering + MCP
#                       Position 21 : MIP-VS + Poisson + Clustering + LASSO + New Tuning Parameter
#                       Position 22 : MIP-VS + Nonlinear + Clustering + LASSO + New Tuning Parameter
#
#   alpha            :  FDR
#   trial.index      :  index for number of simulated replications      
#
#  (METHODOLOGIES)
#   n_folds          :  number of cross-validation 
#   n_subset         :  (MIP & MIP-VS) the number of subsets chosen at random to compute the Min and Max statistics, e.g n_subset=100
#   K.vec            :  a vecotor of number of clusters to be determined, e.g. K.vec <- c(2,3)
#   X.type           :  type of design matrix 
#   output.dir       :  directory for output 
#
# Returns: 
#   decision         :  matrix of decisions  
#   time             :  matrix of time 
#   clean.pos        :  indices of clean observations 

sim.platform.parallel <- function(X, y, family, trial.index, 
                                  detection.switch, 
                                  K.vec, 
                                  X.type, 
                                  n_folds,
                                  n_subset,
                                  alpha,
                                  perturb.type, 
                                  y.infl.pos, 
                                  x.infl.pos, 
                                  y.perturb, 
                                  x.perturb,
                                  output.dir){
  
  #dimension 
  n <- nrow(X)
  p <- ncol(X)
  
  subset.size <- n/2 
  
  ##############################################################################
  # string
  #e.g. "SN1_Poisson2023I_Y5N10_X0N0_IID"
  
  if(all(y.infl.pos==0)==TRUE){
     l1 <- 0 
  }else{
     l1 <- length(y.infl.pos)
  }
 
  if(all(x.infl.pos==0)==TRUE){
     l2 <- 0 
  }else{
     l2 <- length(x.infl.pos)
  }

  string <- paste0("SN", trial.index, "_", perturb.type,"_","Y", l1, "N", y.perturb, 
                   "_","X", l2, "N", x.perturb,"_", X.type) 
  
  ##############################################################################
  
  # names of methods

  names.vec <- c("MIP_VS_LASSO", "MIP_VS_ALASSO", "MIP_VS_SLASSO", "MIP_VS_ENET", "MIP_VS_SCAD", "MIP_VS_MCP",
                 "MIP_VS_Clus_LASSO","MIP_VS_Clus_ALASSO", "MIP_VS_Clus_SLASSO", "MIP_VS_Clus_ENET", "MIP_VS_Clus_SCAD", "MIP_VS_Clus_MCP",
                 "MIP",
                 "DF",
                 "MIPVS(glmLASSO)", "MIPVS(glmSCAD)", "MIPVS(glmMCP)",
                 "ClusMIP(glmLASSO)", "ClusMIP(glmSCAD)", "ClusMIP(glmMCP)",
                 "ClusMIP(PoisLASSO)",
                 "ClusMIP(NonLinearLASSO)")
  
  #extraction 
  method.in.use.pos <- which(detection.switch!=0)
  names.method.in.use <- names.vec[method.in.use.pos]
  
  #decision matrix 
  decision.mat <- matrix(0, nrow=n, ncol=length(method.in.use.pos))
  decision.mat <- as.data.frame(decision.mat)
  colnames(decision.mat) <- names.method.in.use
  
  #time 
  time.mat <- matrix(0, nrow=1, ncol=length(method.in.use.pos))
  time.mat <- as.data.frame(time.mat)
  colnames(time.mat) <- names.method.in.use
  
  #clean data positions 
  clean.pos.mat <- matrix(1, nrow=1, ncol=n)
  
  ##############################################################################
  
  for(b in 1:1){
    
    ############################################################################ 
    
    #full, contaminated dataset 
    X.tilde <- X
    Y.tilde <- y
    Xy.tilde <- as.data.frame(cbind(X.tilde,Y.tilde))
    names(Xy.tilde)<-c(paste("X",1:p,sep=""),"Y")
    
    #clustering
    if(family=="gaussian" | family=="poisson"){
       clean.pos <- cluster.est(df=Xy.tilde, K.vec) 
    }
    
    ############################################################################
    #MIP + VS + LASSO
    #S/N: 1 
    if(detection.switch[1]==1){
       Sys.time()
       tic(paste0("b is ", b, ": MIP + VS + LASSO"))
       tmp.decision.vec <- mult.detection(X=X.tilde, y=Y.tilde, n_subset, subset.size, n_folds, method.switch=c(1,0,0,0,0,0), alpha)$dec.lasso
       exectime <- toc()
       decision.mat$MIP_VS_LASSO <- tmp.decision.vec 
       time.mat$MIP_VS_LASSO <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       Sys.time()
    }
    
    ############################################################################
    #MIP + VS + Adaptive LASSO
    #S/N: 2 
    if(detection.switch[2]==1){
       Sys.time()
       tic(paste0("b is ", b, ": MIP + VS + Adaptive LASSO"))
       tmp.decision.vec <- mult.detection(X=X.tilde, y=Y.tilde, n_subset, subset.size, n_folds, method.switch=c(0,1,0,0,0,0), alpha)$dec.alasso
       exectime <- toc()
       decision.mat$MIP_VS_ALASSO <- tmp.decision.vec 
       time.mat$MIP_VS_ALASSO <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       Sys.time()
    }
    
    ############################################################################
    #MIP + VS + Scaled LASSO
    #S/N: 3 
    if(detection.switch[3]==1){
       Sys.time()
       tic(paste0("b is ", b, ": MIP + VS + Scaled LASSO"))
       tmp.decision.vec <- mult.detection(X=X.tilde, y=Y.tilde, n_subset, subset.size, n_folds, method.switch=c(0,0,1,0,0,0), alpha)$dec.slasso
       exectime <- toc()
       decision.mat$MIP_VS_SLASSO <- tmp.decision.vec 
       time.mat$MIP_VS_SLASSO <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       Sys.time()
    }
    
    ############################################################################
    #MIP + VS + Elastic Net 
    #S/N: 4 
    if(detection.switch[4]==1){
       Sys.time()
       tic(paste0("b is ", b, ": MIP + VS + Elastic Net"))
       tmp.decision.vec <- mult.detection(X=X.tilde, y=Y.tilde, n_subset, subset.size, n_folds, method.switch=c(0,0,0,1,0,0), alpha)$dec.enet
       exectime <- toc()
       decision.mat$MIP_VS_ENET <- tmp.decision.vec 
       time.mat$MIP_VS_ENET <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       Sys.time()
    }
    
    ############################################################################
    #MIP + VS + SCAD
    #S/N: 5 
    if(detection.switch[5]==1){
       Sys.time()
       tic(paste0("b is ", b, ": MIP + VS + SCAD"))
       tmp.decision.vec <- mult.detection(X=X.tilde, y=Y.tilde, n_subset, subset.size, n_folds, method.switch=c(0,0,0,0,1,0), alpha)$dec.vec.scad
       exectime <- toc()
       decision.mat$MIP_VS_SCAD <- tmp.decision.vec 
       time.mat$MIP_VS_SCAD <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       Sys.time()
    }
    
    ############################################################################
    #MIP + VS + MCP
    #S/N: 6 
    if(detection.switch[6]==1){
       Sys.time()
       tic(paste0("b is ", b, ": MIP + VS + MCP"))
       tmp.decision.vec <- mult.detection(X=X.tilde, y=Y.tilde, n_subset, subset.size, n_folds, method.switch=c(0,0,0,0,0,1), alpha)$dec.vec.mcp
       exectime <- toc()
       decision.mat$MIP_VS_MCP <- tmp.decision.vec 
       time.mat$MIP_VS_MCP <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       Sys.time()
    }
    
    ############################################################################
    #MIP + VS + Clustering + LASSO
    #S/N: 7 
    if(detection.switch[7]==1){
       tic(paste0("b is ", b, ": MIP + VS + Clustering + LASSO"))
       tmp.decision.vec <- mult.cluster(X=X.tilde, Y=Y.tilde, n_folds, 
                                        method.switch=c(1,0,0,0,0,0,0,0), 
                                        dat.family=family, clean.pos, alpha)$dec.lasso
       exectime <- toc()
       decision.mat$MIP_VS_Clus_LASSO <- tmp.decision.vec
       time.mat$MIP_VS_Clus_LASSO <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       print(Sys.time())
    }
    
    ############################################################################
    #MIP + VS + Clustering + Adaptive LASSO
    #S/N: 8 
    if(detection.switch[8]==1){
       tic(paste0("b is ", b, ": MIP + VS + Clustering + Adaptive LASSO"))
       tmp.decision.vec <- mult.cluster(X=X.tilde, Y=Y.tilde, n_folds, method.switch=c(0,1,0,0,0,0,0,0), dat.family=family, clean.pos, alpha)$dec.alasso
       exectime <- toc()
       decision.mat$MIP_VS_Clus_ALASSO <- tmp.decision.vec
       time.mat$MIP_VS_Clus_ALASSO <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       print(Sys.time())
    }
    
    ############################################################################
    #MIP + VS + Clustering + Scaled LASSO
    #S/N: 9 
    if(detection.switch[9]==1){
       tic(paste0("b is ", b, ": MIP + VS + Clustering + Scaled LASSO"))
       tmp.decision.vec <- mult.cluster(X=X.tilde, Y=Y.tilde, n_folds, method.switch=c(0,0,1,0,0,0,0,0), dat.family=family, clean.pos, alpha)$dec.slasso
       exectime <- toc()
       decision.mat$MIP_VS_Clus_SLASSO <- tmp.decision.vec 
       time.mat$MIP_VS_Clus_SLASSO <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       print(Sys.time())
    }
    
    ############################################################################
    #MIP + VS + Clustering + Elastic Net 
    #S/N: 10 
    if(detection.switch[10]==1){
       tic(paste0("b is ", b, ": MIP + VS + Clustering + Elastic Net"))
       tmp.decision.vec <- mult.cluster(X=X.tilde, Y=Y.tilde, n_folds, 
                                        method.switch=c(0,0,0,1,0,0,0,0), 
                                        dat.family=family, clean.pos, alpha)$dec.enet
       exectime <- toc()
       decision.mat$MIP_VS_Clus_ENET <- tmp.decision.vec 
       time.mat$MIP_VS_Clus_ENET <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       print(Sys.time())
    }
    
    ############################################################################
    #MIP + VS + Clustering + SCAD 
    #S/N: 11 
    if(detection.switch[11]==1){
       tic(paste0("b is ", b, ": MIP + VS + Clustering + SCAD"))
       tmp.decision.vec <- mult.cluster(X=X.tilde, Y=Y.tilde, n_folds, 
                                        method.switch=c(0,0,0,0,1,0,0,0), 
                                        dat.family=family, clean.pos, alpha)$dec.scad
       exectime <- toc()
       decision.mat$MIP_VS_Clus_SCAD <- tmp.decision.vec 
       time.mat$MIP_VS_Clus_SCAD <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       print(Sys.time())
    }
    
    ############################################################################
    #MIP + VS + Clustering + MCP
    #S/N: 12 
    if(detection.switch[12]==1){
       tic(paste0("b is ", b, ": MIP + VS + Clustering + MCP"))
       tmp.decision.vec <- mult.cluster(X=X.tilde, Y=Y.tilde, n_folds, 
                                        method.switch=c(0,0,0,0,0,1,0,0), 
                                        dat.family=family, clean.pos, alpha)$dec.mcp
       exectime <- toc()
       decision.mat$MIP_VS_Clus_MCP <- tmp.decision.vec
       time.mat$MIP_VS_Clus_MCP <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       print(Sys.time())
    }
  
    ############################################################################
    #MIP
    #S/N: 13  
    if(detection.switch[13]==1){
       tic(paste0("b is ", b, ": MIP"))
       infl.pos.MIP <- MIP(X=as.matrix(X.tilde), Y=t(as.matrix(Y.tilde)), n, p, q=1, n_subset, subset.size, ep = 0.1, alpha)$inf_setfinal
       infl.pos.MIP <- sort(infl.pos.MIP, decreasing=FALSE)
       exectime <- toc()
       decision.mat$MIP[infl.pos.MIP] <- 1  
       time.mat$MIP <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       print(Sys.time())
    }
    
    ############################################################################
    #DF Measure 
    #S/N: 14  
    if(detection.switch[14]==1){
       tic(paste0("b is ", b, ": DF"))
       tmp.decision.vec <- df.measure(X=X.tilde, Y=Y.tilde)
       exectime <- toc()
       decision.mat$DF <- tmp.decision.vec 
       time.mat$DF <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       print(Sys.time())
    }
    
    ############################################################################
    #MIP + VS + Binary + LASSO
    #S/N: 15 
    if(detection.switch[15]==1){
       tic(paste0("b is ", b, ": MIP + VS + Binary + LASSO"))
       tmp.binary.lasso <- mult.detection(X=X.tilde, Y=Y.tilde, n_subset, subset.size, n_folds, method.switch=c(1,0,0,0,0,0), dat.family=family, alpha)
       exectime <- toc()
       #decision.mat$"MIPVS(glmLASSO)" <- tmp.binary.lasso$dec.lasso 
       #myList.unfiltered.decision$"MIPVS(glmLASSO)"[b,] <- tmp.binary.lasso$unfiltered.dec.lasso 
       #time.mat$"MIPVS(glmLASSO)"[b] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       print(Sys.time())
    }
    
    ############################################################################
    #MIP + VS + Binary + SCAD
    #S/N: 16 
    if(detection.switch[16]==1){
       tic(paste0("b is ", b, ": MIP + VS + Binary + SCAD"))
       tmp.binary.scad <- mult.detection(X=X.tilde, Y=Y.tilde, n_subset, subset.size, n_folds, method.switch=c(0,0,0,0,1,0), dat.family=family, alpha)
       exectime <- toc()
       #myList.decision$"MIPVS(glmSCAD)"[b,] <- tmp.binary.scad$dec.scad 
       #myList.unfiltered.decision$"MIPVS(glmSCAD)"[b,] <- tmp.binary.scad$unfiltered.dec.scad 
       #time.mat$"MIPVS(glmSCAD)"[b] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       print(Sys.time())
    }
    
    ############################################################################
    #MIP + VS + Binary + MCP
    #S/N: 17 
    if(detection.switch[17]==1){
       tic(paste0("b is ", b, ": MIP + VS + Binary + MCP"))
       tmp.binary.mcp <- mult.detection(X=X.tilde, Y=Y.tilde, n_subset, subset.size, n_folds, method.switch=c(0,0,0,0,0,1), dat.family=family, alpha)
       exectime <- toc()
       #myList.decision$"MIPVS(glmMCP)"[b,] <- tmp.binary.mcp$dec.mcp
       #myList.unfiltered.decision$"MIPVS(glmMCP)"[b,] <- tmp.binary.mcp$unfiltered.dec.mcp 
       #time.mat$"MIPVS(glmMCP)"[b] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       print(Sys.time())
    }
    
    ############################################################################
    #MIP + VS + Clustering + Binary + LASSO  
    #S/N: 18  
    if(detection.switch[18]==1){
       tic(paste0("b is ", b, ": MIP + VS + Clustering + Binary + LASSO"))
       tmp.decision.vec <- mult.cluster(X=X.tilde, Y=Y.tilde, n_folds, method.switch=c(1,0,0,0,0,0), dat.family=family, clean.pos, alpha)$dec.lasso
       exectime <- toc()
       #myList.decision$"ClusMIP(glmLASSO)"[b,] <- tmp.decision.vec
       #time.mat$"ClusMIP(glmLASSO)"[b] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
       print(Sys.time())
    }
    
    ############################################################################
    #MIP + VS + Clustering + Binary + SCAD  
    #S/N: 19  
    if(detection.switch[19]==1){
      tic(paste0("b is ", b, ": MIP + VS + Clustering + Binary + SCAD"))
      tmp.decision.vec <- mult.cluster(X=X.tilde, Y=Y.tilde, n_folds, method.switch=c(0,0,0,0,1,0), dat.family=family, clean.pos, alpha)$dec.scad
      exectime <- toc()
      #myList.decision$"ClusMIP(glmSCAD)"[b,] <- tmp.decision.vec
      #time.mat$"ClusMIP(glmSCAD)"[b] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
      print(Sys.time())
    }
    
    ############################################################################
    #MIP + VS + Clustering + Binary + MCP  
    #S/N: 20  
    if(detection.switch[20]==1){
      tic(paste0("b is ", b, ": MIP + VS + Clustering + Binary + MCP"))
      tmp.decision.vec <- mult.cluster(X=X.tilde, Y=Y.tilde, n_folds, method.switch=c(0,0,0,0,0,1), dat.family=family, clean.pos, alpha)$dec.mcp
      exectime <- toc()
      myList.decision$"ClusMIP(glmMCP)"[b,] <- tmp.decision.vec
      time.mat$"ClusMIP(glmMCP)"[b] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
      print(Sys.time())
    }
    
    ############################################################################
    #MIP + VS + Poisson + Clustering + LASSO + New Tuning Parameter
    #S/N: 21 
    if(detection.switch[21]==1){
      tic(paste0("b is ", b, ": MIP + VS + Poisson + Clustering + LASSO + New Tuning Parameter"))
      tmp.decision.vec <- mult.cluster(X=X.tilde, Y=Y.tilde, n_folds, 
                                       method.switch=c(0,0,0,0,0,0,1,0), 
                                       dat.family=family, clean.pos, alpha)$dec.lasso.pois
      exectime <- toc()
      decision.mat$"ClusMIP(PoisLASSO)" <- tmp.decision.vec
      time.mat$"ClusMIP(PoisLASSO)" <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
      print(Sys.time())
    }
    
    ############################################################################
    #MIP + VS + Nonlinear + Clustering + LASSO + New Tuning Parameter
    #S/N: 22 
    if(detection.switch[22]==1){
      tic(paste0("b is ", b, ": MIP + VS + Nonlinear + Clustering + LASSO + New Tuning Parameter"))
      tmp.decision.vec <- mult.cluster(X=X.tilde, Y=Y.tilde, n_folds, method.switch=c(0,0,0,0,0,0,0,1), dat.family=family, clean.pos, alpha)$dec.lasso.nlinear
      exectime <- toc()
      myList.decision$"ClusMIP(NonLinearLASSO)"[b,] <- tmp.decision.vec
      time.mat$"ClusMIP(NonLinearLASSO)"[b] <- as.numeric(exectime$toc) - as.numeric(exectime$tic)
      print(Sys.time())
    }
    
  } 
  
  ##############################################################################
  #matrix of decisions 
  write.csv(decision.mat, paste0(output.dir,string,"_dec",".csv"))

  #matrix of time
  write.csv(time.mat, paste0(output.dir,string,"_time",".csv"))
  
  #clean data positions
  write.csv(clean.pos.mat, paste0(output.dir,string,"_cleanPos",".csv"))

  #############################################################################

  returnList <- list("decision"  = decision.mat, 
                     "time"      = time.mat,
                     "clean.pos" = clean.pos.mat)
  
  return(returnList)  
  
}
