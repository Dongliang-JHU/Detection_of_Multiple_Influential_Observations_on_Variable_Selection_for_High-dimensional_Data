# Author: Dongliang Zhang (^)
#
# Supervisors: Professor Masoud Asgharian (*)
#              Professor Martin Lindquist (^)
#              Professpr Mei-Cheng Wang (^)
#
# Current Affiliation: (1) Department of Mathematics and Statistics, McGill University, Montreal, Quebec, Canada (*)
#                      (2) Department of Biostatistics, Johns Hopkins University, Baltimore, Maryland, USA (^)
#
# Title of Project: Detection of Influential Cases on Variable Selection 
#
# R Script Purpose: Step 1 - Multiple Influential Points Detection Algorithm with Clustering
#
# Created on: November 12, 2021 
#
# Modified on        :  April 04, 2022
#                       April 05, 2022
#                       April 07, 2022
#                       April 08, 2022
#                       April 10, 2022
#                       April 13, 2022 
#                       April 14, 2022
#                       April 15, 2022
#                       April 16, 2022 
#                       April 17, 2022
#                       April 27, 2022 
#                       May 10, 2022 
#                       May 15, 2022 
#                       May 16, 2022 
#                       May 20, 2022 
#                       June 02, 2022 
#                       August 25, 2022 
#
# Dependent R Scripts: lasso_fit.r, 
#                      enet_fit.r, 
#                      normalDecision.r; 
#
# R Packages         : glmnet, 
#                      scalreg,
#                      ncvreg,
#                      SIS,
#                      mpath; 
#                      
#############################################################################################################################
# Function: mult.cluster
#
# Args:
#
#   X                :  potentially contaminated design matrx 
#   Y                :  potentially contaminated response vector  
#   n_folds          :  nfold cross-validation for LASSO and elastic net 
#   method.switch    :  a switch control (0 or 1) for model selection methods 
#                       Position 1: LASSO
#                       Position 2: adaptive LASSO 
#                       Position 3: scaled LASSO
#                       Position 4: elastic net 
#                       Position 5: SCAD
#                       Position 6: MCP
#                       Position 7: LASSO + Poisson + New Tuning Parameter 
#                       Position 8: LASSO + Nonlinear + New Tuning Parameter 
#   dat.family       :  gaussian or binomial or poisson 
#   clean.pos        :  indices of estimated clean dataset from clustering  
#   alpha            :  FDR 
#   verbose          :  TRUE or FALSE 
#
# Returns: 
#
#   dec.lasso         :  a decision vector (0 or 1) for LASSO
#   dec.alasso        :  a decision vector (0 or 1) for adaptive LASSO
#   dec.slasso        :  a decision vector (0 or 1) for scaled LASSO 
#   dec.enet          :  a decision vector (0 or 1) for elastic net 
#   dec.scad          :  a decision vector (0 or 1) for SCAD 
#   dec.mcp           :  a decision vector (0 or 1) for MCP 
#
# Dependent R functions: mult.filtering  

mult.cluster <- function(X, Y, 
                         n_folds, 
                         method.switch,
                         dat.family,
                         clean.pos,
                         alpha){
  
  n <- nrow(X)
  p <- ncol(X)
  
  #decision vector 
  dec.vec.lasso          <- rep(0, n)
  dec.vec.alasso         <- rep(0, n)
  dec.vec.slasso         <- rep(0, n)
  dec.vec.enet           <- rep(0, n)
  dec.vec.scad           <- rep(0, n)
  dec.vec.mcp            <- rep(0, n) 
  dec.vec.lasso.pois     <- rep(0, n) 
  dec.vec.lasso.nlinear  <- rep(0, n)

  infl.prelim.pos <- setdiff(1:n, clean.pos)
  
  #print(infl.prelim.pos)
  
  ####################################################################################################################
  #LASSO
  if(method.switch[1]==1){
  
     S.lasso <- c()
     
     for(i in 1:length(infl.prelim.pos)){
       
         X.sub <- rbind(X[infl.prelim.pos[i],], X[clean.pos,])
         Y.sub <- c(Y[infl.prelim.pos[i]], Y[clean.pos])  
        
         glm.obj <- cv.glmnet(X.sub, Y.sub, alpha=1, family=dat.family, nfolds=n_folds)   
         coeff <- as.numeric(coef(glm.obj, s="lambda.min")[-1]) 
        
         boo.coeff <- which(coeff==0) 
         boo.org <- rep(0,p) 
         boo.org[boo.coeff] <- 1  
        
         tau.vec <- sapply(1:nrow(X.sub), function(j){
           Y.sub.red <- Y.sub[-j] 
           X.sub.red <- X.sub[-j,] 
           glm.obj.red <- cv.glmnet(X.sub.red, Y.sub.red, alpha=1, family=dat.family, nfolds=n_folds)  
           coeff.red <- as.numeric(coef(glm.obj.red, s="lambda.min")[-1])  
           boo.coeff.red <- which(coeff.red==0) 
           boo.red <- rep(0,p) 
           boo.red[boo.coeff.red] <- 1  
           sum((boo.org-boo.red)^2)})
  
         if(normal.Decision(tau.vec, scaling.factor=0.2, alpha)[1] == 1){
            S.lasso <- c(S.lasso, infl.prelim.pos[i])
         }
        
     }#end of for loop
    
     if(length(S.lasso) > 1){
        dec.vec.lasso[S.lasso] <- 1 
     }
    
  }#end of LASSO 
  
  ####################################################################################################################
  #Adaptive LASSO
  if(method.switch[2]==1){

    S.alasso <- c()
    
    for(i in 1:length(infl.prelim.pos)){
      
      X.sub <- rbind(X[infl.prelim.pos[i],], X[clean.pos,])
      Y.sub <- c(Y[infl.prelim.pos[i]], Y[clean.pos])  
    
      ridge.cv    <- cv.glmnet(x=X.sub, y=Y.sub, type.measure="mse", nfolds=n_folds, alpha=0)
      ridge.coef  <- as.numeric(coef(ridge.cv, s=ridge.cv$lambda.min))[-1]
      
      alasso.cv   <- cv.glmnet(x=X.sub, y=Y.sub, type.measure="mse", nfolds=n_folds, alpha=1, penalty.factor=1/abs(ridge.coef), keep=TRUE)
      alasso.coef <- as.numeric(coef(alasso.cv, s=alasso.cv$lambda.min))[-1]
      
      boo.coeff   <- which(alasso.coef==0) 
      boo.org     <- rep(0, p) 
      boo.org[boo.coeff] <- 1  
      
      tau.vec <- sapply(1:nrow(X.sub), function(j){
        Y.sub.red <- Y.sub[-j] 
        X.sub.red <- X.sub[-j,] 
        
        ridge.cv.red <- cv.glmnet(x=X.sub.red, y=Y.sub.red, type.measure="mse", nfolds=n_folds, alpha=0)
        ridge.coef.red <- as.numeric(coef(ridge.cv.red, s=ridge.cv.red$lambda.min))[-1]
        
        alasso.cv.red <- cv.glmnet(x=X.sub.red, y=Y.sub.red, type.measure="mse", nfolds=n_folds, alpha=1, penalty.factor=1/abs(ridge.coef.red), keep=TRUE)
        alasso.coef.red <- as.numeric(coef(alasso.cv.red, s=alasso.cv.red$lambda.min))[-1]
        
        boo.coeff.red <- which(alasso.coef.red==0) 
        boo.red <- rep(0, p) 
        boo.red[boo.coeff.red] <- 1  
        sum((boo.org-boo.red)^2)})
      
      if(normal.Decision(tau.vec, scaling.factor=0.2, alpha)[1] == 1){
        S.alasso <- c(S.alasso, infl.prelim.pos[i])
      }
      
    }#end of for loop
    
    if(length(S.alasso) > 1){
      dec.vec.alasso[S.alasso] <- 1 
    }
    
  }#end of Adaptive LASSO 
  
  ####################################################################################################################
  #Scaled LASSO
  if(method.switch[3]==1){
    
     S.slasso <- c() 
    
     for(i in 1:length(infl.prelim.pos)){
    
        X.sub <- rbind(X[infl.prelim.pos[i],], X[clean.pos,])
        Y.sub <- c(Y[infl.prelim.pos[i]], Y[clean.pos])  
        
        slasso.obj <- scalreg(X.sub, Y.sub, lam0 = NULL, LSE = FALSE)
        coeff <- as.numeric(slasso.obj$coefficients)  
        boo.coeff <- which(coeff==0) 
        boo.org <- rep(0,p) 
        boo.org[boo.coeff] <- 1    
        
        tau.vec <- sapply(1:nrow(X.sub), function(j){
          Y.sub.red <- Y.sub[-j]  
          X.sub.red <- X.sub[-j,]   
          slasso.obj.red <- scalreg(X.sub.red, Y.sub.red, lam0 = NULL, LSE = FALSE) 
          coeff.red <- as.numeric(slasso.obj.red$coefficients)  
          boo.coeff.red <- which(coeff.red==0) 
          boo.red <- rep(0,p) 
          boo.red[boo.coeff.red] <- 1  
          sum((boo.org-boo.red)^2) 	  
        })
        
        if(normal.Decision(tau.vec, scaling.factor=0.2, alpha)[1] == 1){
           S.slasso <- c(S.slasso, infl.prelim.pos[i])
        }
        
    }#end of "for" loop 
    
    if(length(S.slasso) > 1){
      dec.vec.slasso[S.slasso] <- 1 
    }
    
  }#end of Scaled LASSO 
  
  ####################################################################################################################
  #Elastic Net 
  if(method.switch[4]==1){
    
     S.enet <- c()  
     
     for(i in 1:length(infl.prelim.pos)){
        
         X.sub <- rbind(X[infl.prelim.pos[i],], X[clean.pos,])
         Y.sub <- c(Y[infl.prelim.pos[i]], Y[clean.pos])  
        
         coeff <- enet.fit(X=X.sub, Y=Y.sub, n_folds, dat.family)$enet.coef
         boo.coeff <- which(coeff==0) 
         boo.org <- rep(0,p) 
         boo.org[boo.coeff] <- 1       
         
         tau.vec <- sapply(1:nrow(X.sub), function(j){
           Y.sub.red <- Y.sub[-j]  
           X.sub.red <- X.sub[-j,]   
           coeff.red <- enet.fit(X=X.sub.red, Y=Y.sub.red, n_folds, dat.family)$enet.coef
           boo.coeff.red <- which(coeff.red==0) 
           boo.red <- rep(0,p) 
           boo.red[boo.coeff.red] <- 1  
          sum((boo.org-boo.red)^2)  
         })
         
         if(normal.Decision(tau.vec, scaling.factor=0.2, alpha)[1]==1){
            S.enet <- c(S.enet, infl.prelim.pos[i])
         }
        
     }#end of for loop
    
     if(length(S.enet) > 1){
        dec.vec.enet[S.enet] <- 1 
     }
    
  }#end of Elastic Net 
  
  ####################################################################################################################  
  #SCAD 
  if(method.switch[5]==1){
    
     S.scad <- c() 
    
     for(i in 1:length(infl.prelim.pos)){
        
         X.sub <- rbind(X[infl.prelim.pos[i],], X[clean.pos,])
         Y.sub <- c(Y[infl.prelim.pos[i]], Y[clean.pos])  
        
         cv.scad <- cv.ncvreg(X.sub, Y.sub, family=dat.family, penalty="SCAD") 
         scad.obj <- cv.scad$fit 
         coeff <- as.numeric(scad.obj$beta[,cv.scad$min][-1])
         boo.coeff <- which(coeff==0) 
         boo.org <- rep(0,p) 
         boo.org[boo.coeff] <- 1      
         
         tau.vec <- sapply(1:nrow(X.sub), function(j){
           Y.sub.red <- Y.sub[-j]  
           X.sub.red <- X.sub[-j,]   
           cv.scad.red <- cv.ncvreg(X.sub.red, Y.sub.red, family=dat.family, penalty="SCAD")  
           scad.obj.red <- cv.scad.red$fit 
           coeff.red <- as.numeric(scad.obj.red$beta[,cv.scad.red$min][-1])
           boo.coeff.red <- which(coeff.red==0) 
           boo.red <- rep(0,p) 
           boo.red[boo.coeff.red] <- 1  
           sum((boo.org-boo.red)^2) 	  
         })
        
         if(normal.Decision(tau.vec, scaling.factor=0.2, alpha)[1] == 1){
            S.scad <- c(S.scad, infl.prelim.pos[i])
         }
        
     }#end of for loop  
    
     if(length(S.scad) > 1){
        dec.vec.scad[S.scad] <- 1 
     }
    
  }#end of SCAD 
  
  ####################################################################################################################  
  #MCP 
  if(method.switch[6]==1){

    S.mcp <- c() 
      
    for(i in 1:length(infl.prelim.pos)){
      
        tau.vec <- vector( mode="numeric", length=(length(clean.pos)+1) ) 
      
        skip_to_next_outer <- FALSE
        
        X.sub <- rbind(X[infl.prelim.pos[i],], X[clean.pos,])
        Y.sub <- c(Y[infl.prelim.pos[i]], Y[clean.pos])  
        
        tryCatch({cv.mcp <- cv.ncvreg(X.sub, Y.sub, family=dat.family, penalty="MCP") 
                  mcp.obj <- cv.mcp$fit},
                  error = function(e){
                  skip_to_next_outer <<- TRUE
                  message(paste0("An error occurred for potential influential observation indexed by ", infl.prelim.pos[i], ":\n"), e)
                  }
        )
        
        if(skip_to_next_outer){
          S.mcp <- c(S.mcp, infl.prelim.pos[i])
          next
        }
      
        coeff <- as.numeric(mcp.obj$beta[,cv.mcp$min][-1])
        
        boo.coeff <- which(coeff==0) 
        boo.org <- rep(0,p) 
        boo.org[boo.coeff] <- 1       
        
        for(j in 1:nrow(X.sub)){
          
            skip_to_next <- FALSE
            
            Y.sub.red <- Y.sub[-j]  
            X.sub.red <- X.sub[-j,]   
            
            tryCatch({cv.mcp.red <- cv.ncvreg(X.sub.red, Y.sub.red, family=dat.family, penalty="MCP")}, 
                      error = function(e){
                      skip_to_next <<- TRUE
                      message(paste0("An error occurred for potential influential observation indexed by ", infl.prelim.pos[i], " when j is ", j, ":\n"), e)
                      }
            )
            
            if(skip_to_next){
              tau.vec[j] <- 0
              next
            }
             
            mcp.obj.red <- cv.mcp.red$fit 
          
            coeff.red <- as.numeric(mcp.obj.red$beta[,cv.mcp.red$min][-1])
            
            boo.coeff.red <- which(coeff.red==0) 
            
            boo.red <- rep(0,p) 
            
            boo.red[boo.coeff.red] <- 1  
            
            tau.vec[j] <- sum((boo.org-boo.red)^2)  	
            
         } 
        
         if(normal.Decision(tau.vec, scaling.factor=0.2, alpha)[1] == 1){
            S.mcp <- c(S.mcp, infl.prelim.pos[i])
         }
      
     }#end of for loop 
    
     if(length(S.mcp) > 1){
        dec.vec.mcp[S.mcp] <- 1 
     } 
    
  }#end of MCP 
  
  ####################################################################################################################
  #LASSO + Poisson + New Tuning Parameter
  if(method.switch[7]==1){
     
     S.lasso.pois <- c()
    
     for(i in 1:length(infl.prelim.pos)){
      
         tau.vec <- vector( mode="numeric", length=(length(clean.pos)+1) )   
         
         skip_to_next_outer <- FALSE
      
         X.sub <- rbind(X[infl.prelim.pos[i],], X[clean.pos,])
         Y.sub <- c(Y[infl.prelim.pos[i]], Y[clean.pos])  
         n.sub <- length(Y.sub)
      
         lam     <- 1.1 * (1/(sqrt(n.sub))) * qnorm(1-alpha/(2*p), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
         lam.red <- 1.1 * (1/(sqrt(n.sub-1))) * qnorm(1-alpha/(2*p), mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
         
         tryCatch({fit.pois <- glmreg(x=X.sub, y=Y.sub, family = "poisson", lambda=lam)},
                  error = function(e){
                  skip_to_next_outer <<- TRUE
                  message(paste0("An error occurred for potential influential observation indexed by ", infl.prelim.pos[i], ":\n"), e)
                  }
         )
         
         if(skip_to_next_outer){
           S.lasso.pois <- c(S.lasso.pois, infl.prelim.pos[i])
           next
         }
         
         coeff <- as.numeric(coef(fit.pois))[-1]
         boo.coeff <- which(coeff==0) 
         boo.org <- rep(0,p) 
         boo.org[boo.coeff] <- 1  
      
         for(j in 1:nrow(X.sub)){
           
           skip_to_next <- FALSE
           
           tryCatch({
             Y.sub.red <- Y.sub[-j] 
             X.sub.red <- X.sub[-j,]
             fit.pois.red <- glmreg(x=X.sub.red, y=Y.sub.red, family = "poisson", lambda=lam.red)
             coeff.red <- as.numeric(coef(fit.pois.red))[-1]
             boo.coeff.red <- which(coeff.red==0) 
             boo.red <- rep(0,p) 
             boo.red[boo.coeff.red] <- 1  
             tau.vec[j] <- sum((boo.org-boo.red)^2)},
             
             error = function(e){
               skip_to_next <<- TRUE
               message(paste0("An error occurred for potential influential observation indexed by ", infl.prelim.pos[i], " when j is ", j, ":\n"), e)
             })
             
             if(skip_to_next){
               tau.vec[j] <- 0
               next
             }
          
         }
      
         if(normal.Decision(tau.vec, scaling.factor=0.2, alpha)[1] == 1){
            S.lasso.pois <- c(S.lasso.pois, infl.prelim.pos[i])
         }
      
     }#end of for loop
    
     if(length(S.lasso.pois) > 1){
        dec.vec.lasso.pois[S.lasso.pois] <- 1 
     }
    
  }#end of LASSO + Poisson + New Tuning Parameter 
  
  ####################################################################################################################
  #LASSO + Nonlinear + New Tuning Parameter
  if(method.switch[8]==1){
    
     S.lasso.nlinear <- c()
    
     for(i in 1:length(infl.prelim.pos)){
      
         X.sub <- rbind(X[infl.prelim.pos[i],], X[clean.pos,])
         Y.sub <- c(Y[infl.prelim.pos[i]], Y[clean.pos])  
        
         lambda.vec <- seq(from=n^0.4, to=n^0.9, length.out=100)
         lambda.vec.red <- seq(from=(n-1)^0.4, to=(n-1)^0.9, length.out=100)
         
         glm.obj <- cv.glmnet(x=X.sub, y=Y.sub, alpha=1, family=dat.family, nfolds=n_folds, lambda=lambda.vec)
         coeff <- as.numeric(coef(glm.obj, s="lambda.min")[-1]) 
         
         boo.coeff <- which(coeff==0) 
         boo.org <- rep(0,p) 
         boo.org[boo.coeff] <- 1  
      
         tau.vec <- sapply(1:nrow(X.sub), function(j){
           Y.sub.red <- Y.sub[-j] 
           X.sub.red <- X.sub[-j,]
           glm.obj.red <- cv.glmnet(x=X.sub.red, y=Y.sub.red, alpha=1, family=dat.family, nfolds=n_folds, lambda=lambda.vec.red)
           coeff.red <- as.numeric(coef(glm.obj.red, s="lambda.min")[-1]) 
           boo.coeff.red <- which(coeff.red==0) 
           boo.red <- rep(0,p) 
           boo.red[boo.coeff.red] <- 1  
           sum((boo.org-boo.red)^2) 
         })
      
         if(normal.Decision(tau.vec, scaling.factor=0.2, alpha)[1] == 1){
            S.lasso.nlinear <- c(S.lasso.nlinear, infl.prelim.pos[i])
         }
      
     }#end of for loop
    
     if(length(S.lasso.nlinear) > 1){
        dec.vec.lasso.nlinear[S.lasso.nlinear] <- 1 
     }
    
  }#end of LASSO + Nonlinear + New Tuning Parameter 
  
  #################################################################################################################### 
  
  returnList <- list("dec.lasso"         = dec.vec.lasso,
                     "dec.alasso"        = dec.vec.alasso,
                     "dec.slasso"        = dec.vec.slasso,
                     "dec.enet"          = dec.vec.enet, 
                     "dec.scad"          = dec.vec.scad, 
                     "dec.mcp"           = dec.vec.mcp,
                     "dec.lasso.pois"    = dec.vec.lasso.pois,
                     "dec.lasso.nlinear" = dec.vec.lasso.nlinear)
  
  return(returnList)  
  
}
