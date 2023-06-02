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
# R Script Purpose: Simulation Plot of (1) Power and FPR, and (2) Time   
# 
# Created on  : March 22, 2023 
#
# Modified on : March 22, 2023 
#
# Dependent R Scripts:  
#
# R Packages : ggplot2,
#              stringr, 
#              gridExtra, 
#              lemon (grid_arrange_shared_legend)  
#              stringr
#
##########################################################################################################################################
# Function: sim.plot 
#
# Args:
# 
#   df.list          :  a list of dataframes containing power and FPR  
#   infl.prop        :  a vector of proportion of influential observations 
#   perturb.type     :  perturbation type
#   y.infl.pos       :  indices of influential response 
#   x.infl.pos       :  indices of influential explanatory variables 
#   y.perturb        :  magnitude of perturbation on response vector  
#   x.perturb        :  magnitude of perturbation on explanatory variables  
#
# Returns: 
#   power.dat        :  dataframe of power 
#   FPR.dat          :  dataframe of FPR 
#   time.dat         :  dataframe of time 

sim.plot <- function(source.dir,
                     output.dir,
                     X.type, 
                     perturb.type, 
                     y.perturb, 
                     x.perturb){
  
  #e.g. "Zhao2015I_Y10_X0_IID"
  myString <- paste0(perturb.type,"_","Y",y.perturb, "_","X", x.perturb,"_", X.type) 
  
  #list.l <- length(df.list)
  
  #combined dataframe
  #if(list.l==1){
  #  df <- t(data.frame(Reduce(rbind, df.list)))
  #}else{
  #  df <- data.frame(Reduce(rbind, df.list))
  #}
  
  #names.str <- as.character(df[,1]) 
  
  #names.str.new <- str_replace_all(names.str, c("MIP_VS_LASSO"            = "MIP-VS-LASSO", 
  #                                              "MIP_VS_ALASSO"           = "MIP-VS-ALASSO",
  #                                              "MIP_VS_SASSO"            = "MIP-VS-SLASSO", 
  #                                              "MIP_VS_ENET"             = "MIP-VS-ENET",
  #                                              "MIP_VS_SCAD"             = "MIP-VS-SCAD", 
  #                                              "MIP_VS_MCP"              = "MIP-VS-MCP",
  #                                              "MIP_VS_Clus_LASSO"       = "ClusMIP(LASSO)", 
  #                                              "MIP_VS_Clus_ALASSO"      = "ClusMIP(ALASSO)",
  #                                              "MIP_VS_Clus_SLASSO"      = "ClusMIP(SLASSO)", 
  #                                              "MIP_VS_Clus_ENET"        = "ClusMIP(ENET)",
  #                                              "MIP_VS_Clus_SCAD"        = "ClusMIP(SCAD)", 
  #                                              "MIP_VS_Clus_MCP"         = "ClusMIP(MCP)",
  #                                              "ClusMIP\\(PoisLASSO\\)"  = "ClusMIP\\(tunedLASSO\\)"
  #))
  
  #df[,1] <- as.factor(names.str.new)
  
  #method.names <- as.character(df[,1])[1:nrow(df.list[[1]])]
  
  ##########################################################################################
  #dataframe of power  
  #df.power <- data.frame(Methods=rep(method.names,list.l),
  #                       Power=df$Power,
  #                       Proportion=rep(infl.prop, each=nrow(df.list[[1]])))
  #df.power$Methods <- factor(df.power$Methods, levels = method.names)
  #df.power$Proportion <- factor(df.power$Proportion, levels = infl.prop)
  
  ##########################################################################################
  #dataframe of FPR 
  #df.fpr <- data.frame(Methods=rep(method.names,list.l),
  #                     FPR=df$FPR,
  #                     Proportion=rep(infl.prop, each=nrow(df.list[[1]])))
  #df.fpr$Methods <- factor(df.fpr$Methods, levels = method.names)
  #df.fpr$Proportion <- factor(df.fpr$Proportion, levels = infl.prop)
  
  ##########################################################################################
  #dataframe of time
  #df.time <- data.frame(Methods=rep(method.names,list.l),
  #                      Time=df$Time,
  #                      Proportion=rep(infl.prop, each=nrow(df.list[[1]])))
  #df.time$Methods <- factor(df.time$Methods, levels = method.names)
  #df.time$Proportion <- factor(df.time$Proportion, levels = infl.prop)
  
  ##########################################################################################
  #import  
  
  df.power <- read.csv(paste0(source.dir,"Power_",myString,".csv"))
  df.fpr <- read.csv(paste0(source.dir,"FPR_",myString,".csv"))
  df.time <- read.csv(paste0(source.dir,"Time_",myString,".csv"))
 
  names.str <- as.character(df.power$Methods) 
  names.str.new <- str_replace_all(names.str, c("DF"                    = "DF\\(LASSO\\)",
                                                "ClusMIP\\(LASSO\\)"    = "ClusMIP\\(LASSO\\)",
                                                "ClusMIP\\(SLASSO\\)"   = "ClusMIP\\(SLASSO\\)",
                                                "ClusMIP\\(ENET\\)"     = "ClusMIP\\(ENET\\)",
                                                "ClusMIP\\(SCAD\\)"     = "ClusMIP\\(SCAD\\)",
                                                "ClusMIP\\(MCP\\)"      = "ClusMIP\\(MCP\\)",
                                                "ClusMIP\\(tunedLASSO\\)"   = "ClusMIP(WLASSO)")) 
  
  df.power$Methods <- as.factor(names.str.new)
  method.names <- as.character(unique(df.power$Methods))
  df.power$Methods <- factor(df.power$Methods, levels = method.names)

  df.fpr$Methods <- as.factor(names.str.new)
  df.fpr$Methods <- factor(df.fpr$Methods, levels = method.names)
  
  df.time$Methods <- as.factor(names.str.new)
  df.time$Methods <- factor(df.time$Methods, levels = method.names)

  ##########################################################################################
  infl.prop <- c("5%","10%","15%","20%")
  
  size.tmp <- 32
    
  jpeg(filename=paste0(output.dir, myString, ".jpeg"), width=2200, height=1000) 
  
  fig1 <- ggplot(df.power, aes(x=Methods, y=Power, group=Proportion)) +
    geom_line(linetype = "dashed") + 
    #geom_point(aes(shape=Methods, color=Proportion), size = 4)  +
    geom_point(aes(color=Proportion), size = 4)  +
    xlab("") + ylab("Power") + 
    theme(axis.title.x = element_text(color="black", size=size.tmp, face="bold"),
          axis.title.y = element_text(color="black", size=size.tmp, face="bold"),
          axis.text.x = element_text(color="black",  size=size.tmp, angle=45, vjust=0.5),
          axis.text.y = element_text(color="black",  size=size.tmp),
          legend.title = element_text(size=size.tmp, face="bold"),
          legend.text = element_text(size=size.tmp),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.margin = unit(c(3,3,3,3), "lines")) + 
    scale_y_continuous(breaks = seq(0, 1.0, 0.1)) + 
    coord_cartesian(ylim = c(0, 1.0)) +   
    scale_colour_discrete(name="Contamination Proportion", labels=infl.prop) #+ 
    #scale_shape_discrete(name="Methods", labels=method.names)
  
  fig2 <- ggplot(df.fpr, aes(x=Methods, y=FPR, group=Proportion)) +
    geom_line(linetype = "dashed") + 
    #geom_point(aes(shape=Methods, color=Proportion), size = 4)  +
    geom_point(aes(color=Proportion), size = 4)  +
    xlab("") + ylab("FPR") + 
    theme(axis.title.x = element_text(color="black", size=size.tmp, face="bold"),
          axis.title.y = element_text(color="black", size=size.tmp, face="bold"),
          axis.text.x = element_text(color="black", size=size.tmp, angle=45, vjust=0.5),
          axis.text.y = element_text(color="black", size=size.tmp),
          legend.title = element_text(size=size.tmp, face="bold"),
          legend.text = element_text(size=size.tmp),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.margin = unit(c(3,3,3,3), "lines")) + 
    scale_y_continuous(breaks = seq(0, 1.0, 0.1)) + 
    coord_cartesian(ylim = c(0, 1.0)) +   
    scale_colour_discrete(name="Contamination Proportion", labels=infl.prop) #+ 
    #scale_shape_discrete(name="Methods", labels=method.names)
  
  fig3 <- ggplot(df.time, aes(x=Methods, y=Time, group=Proportion)) +
    geom_line(linetype = "dashed") + 
    #geom_point(aes(shape=Methods, color=Proportion), size = 4)  +
    geom_point(aes(color=Proportion), size = 4)  +
    xlab("") + ylab("Time in Seconds") + 
    theme(axis.title.x = element_text(color="black", size=size.tmp, face="bold"),
          axis.title.y = element_text(color="black", size=size.tmp, face="bold"),
          axis.text.x = element_text(color="black", size=size.tmp, angle=45, vjust=0.5),
          axis.text.y = element_text(color="black", size=size.tmp),
          legend.title = element_text(size=size.tmp, face="bold"),
          legend.text = element_text(size=size.tmp),
          panel.border = element_rect(color = "black", fill = NA, size = 1),
          plot.margin = unit(c(3,3,3,3), "lines")) + 
    #scale_y_continuous(breaks = seq(0, 1.0, 0.1)) + 
    #coord_cartesian(ylim = c(0, 1.0)) +   
    scale_colour_discrete(name="Contamination Proportion", labels=infl.prop) #+ 
    #scale_shape_discrete(name="Methods", labels=method.names)
  
  final <- ggarrange(fig1, fig2, fig3,  
                     #labels = c("(a): low-dimensional", 
                     #            "(b): high-dimensional"),
                     ncol = 3, nrow = 1, vjust=-0.05, hjust=-0.60, common.legend = TRUE, font.label=list(color="black",                       size=42, face="bold"))
  
  annotate_figure(final, top = text_grob("", color = "black", face = "bold", size = 42))
  
  print(final)
  
  dev.off()
  
  
  ##########################################################################################
  
  #returnList <- list("power.dat"  = df.power,
  #                   "FPR.dat"    = df.fpr, 
  #                   "time.dat"   = df.time)
  
}