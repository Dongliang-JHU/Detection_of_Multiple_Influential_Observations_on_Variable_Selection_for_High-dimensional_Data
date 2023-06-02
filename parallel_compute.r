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
# R Script Purpose: Parallel Computing of Power, FPR and Time 
# 
# Created on  : May 31, 2022 
#
# Modified on : June 01, 2022 
#               June 04, 2022 
#
# Dependent R Scripts:  testsEval.r; 
#
# R Packages :          rlist (saving a list);


parallel.computing <- function(n, n_sim, detection.switch, 
                               perturb.type, 
                               y.infl.pos, y.perturb, x.infl.pos, x.perturb, 
                               X.type){
  
  #########################################################################################################################################
  names.vec <- c("MIP_VS_LASSO", "MIP_VS_ALASSO", "MIP_VS_SLASSO", "MIP_VS_ENET", "MIP_VS_SCAD", "MIP_VS_MCP",
                 "MIP_VS_Clus_LASSO","MIP_VS_Clus_ALASSO", "MIP_VS_Clus_SLASSO", "MIP_VS_Clus_ENET", "MIP_VS_Clus_SCAD", "MIP_VS_Clus_MCP",
                 "MIP",
                 "DF",
                 "MIPVS(glmLASSO)", "MIPVS(glmSCAD)", "MIPVS(glmMCP)",
                 "ClusMIP(glmLASSO)", "ClusMIP(glmSCAD)", "ClusMIP(glmMCP)",
                 "ClusMIP(PoisLASSO)",
                 "ClusMIP(NonLinearLASSO)")
  
  method.in.use.pos <- which(detection.switch!=0)
  
  names.method.in.use <- names.vec[method.in.use.pos]
  
  length.method <- length(which(detection.switch!=0))
  
  #string
  #e.g. "Poisson2023I_Y5%10_X0%0_IID"
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
  
  #string <- paste0(perturb.type,"_","Y",format(round((l1/n*100),2),nsmall=2),"%",y.perturb, "_","X", round(l2/n*100), "%", x.perturb,"_",X.type)
  
  string <- paste0(perturb.type,"_","Y",round(l1/n*100),"%",y.perturb, "_","X", round(l2/n*100), "%", x.perturb,"_",X.type)
  
  #########################################################################################################################################
  dec.list <- vector(mode="list", length=n_sim)
  
  time.mat <- matrix(nrow=n_sim, ncol=length.method)
  
  Y.mat <- matrix(nrow=n_sim, ncol=n)
  
  #########################################################################################################################################
  #import 
  for(i in 1:n_sim){
    #dec.list[[i]] <- as.matrix(read.csv(paste0("SN", i, "_", perturb.type,"_","Y", l1, "N", y.perturb, 
    #                                 "_","X", l2, "N", x.perturb,"_", X.type,"_dec.csv")))[,-1]
    #time.mat[i,] <- as.matrix(read.csv(paste0("SN", i, "_", perturb.type,"_","Y", l1, "N", y.perturb, 
    #                                "_","X", l2, "N", x.perturb,"_", X.type,"_time.csv")))[,-1]
    dec.list[[i]] <- as.matrix(read.csv(paste0("SN", i, "_", perturb.type,"_","Y",(l1/n*100),"%",y.perturb, 
                                               "_","X", round(l2/n*100), "%", x.perturb,"_",X.type,"_dec.csv")))[,-1]
    time.mat[i,] <- as.matrix(read.csv(paste0("SN", i, "_", perturb.type,"_","Y",(l1/n*100),"%",y.perturb, 
                                              "_","X", round(l2/n*100), "%", x.perturb,"_",X.type,"_time.csv")))[,-1]
    Y.mat[i,] <- as.matrix(read.csv(paste0("SN", i, "_", perturb.type,"_","Y",(l1/n*100),"%",y.perturb, 
                                           "_","X", round(l2/n*100), "%", x.perturb,"_",X.type,"_Y.csv")))[,-1]
  }
  
  #decision list formation 
  myList.decision <- rep(list(matrix(0, nrow=n_sim, ncol=n)), length(names.method.in.use))
  names(myList.decision) <- names.method.in.use
  
  for(j in 1:length(names.method.in.use)){
    tmp <- mapply('[', dec.list, TRUE, rep(j,100))
    tmp.mat <- t(as.matrix(tmp))
    colnames(tmp.mat) <- NULL
    myList.decision[[j]] <- tmp.mat
  }
  
  
  #########################################################################################################################################
  #time computation
  time <- as.numeric(colMeans(time.mat))
  time.mat <- data.frame(Methods=names.method.in.use, Time=time)
  
  #power, FPR and time computation  
  list.result <- lapply(myList.decision, testsEval, infl.index=y.infl.pos, noninfl.index=setdiff(1:n,y.infl.pos))
  
  if(length(method.in.use.pos)==1){
    dat.result <- t(data.frame(Reduce(rbind, list.result)))
  }else{
    dat.result <- data.frame(Reduce(rbind, list.result))
  }
  
  dat.result <- data.frame(dat.result, time)
  colnames(dat.result) <- c("Power", "FPR", "FDR", "Time")
  rownames(dat.result) <- names.method.in.use
  dat.result <- dat.result[,-3]
  
  #########################################################################################################################################
  #export 
  #Y
  write.csv(Y.mat, paste0(string,"_Y",".csv"))
  
  #list of decisions 
  list.save(myList.decision, paste0(string,"_dec",".rdata"))
  #C <- list.load(paste0(string,"_dec",".rdata"))
  
  #time
  write.csv(time.mat, paste0(string,"_time",".csv"))
  
  #final results 
  write.csv(dat.result, paste0(string,"_rslt",".csv"))
  
}