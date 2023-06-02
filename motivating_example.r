# Author: Dongliang Zhang (^)
#
# Supervisors: Professor Masoud Asgharian (*)
#              Professor Martin A. Lindquist (^)
#
# Current Affiliation: (1) Department of Mathematics and Statistics, McGill University, Montreal, Quebec, Canada (*)
#                      (2) Department of Biostatistics, Johns Hopkins University, Baltimore, Maryland, USA (^)
#
# Title of Project: Detection of Influential Cases on Variable Selection 
#
# R Script Purpose: Motivating Example 
#
# Created on      : March 20, 2023
#
# Modified on     : March 20, 2023
#
# Dependent R Scripts: dataGenerate.r;  
#
# R Packages         : MASS; 

################################################################################

setwd("/Users/dongliangzhang/Library/CloudStorage/Dropbox/Research/Detection of Influential Observations on Variable Selection/Dongliang Zhang_Influential Diagnostics/Code/Motivating Example/Final")

library(mvtnorm)
library(glmnet)
library(ncvreg)
library(scalreg)
library(rje) #is.subset 
library(ggplot2)
library(ggpubr)
library(bestglm)
library(MASS)
library(rlist)
library(gridExtra)
library(gtable)
library(grid)

source("dataGenerate.r")
source("Simulation I.r")

################################################################################

# low-dimensional
# perturbation model I 
df.low.Zhao2015I <- as.data.frame(read.csv("Example I_Simulation I_all_dat_Zhao2015I.csv")[,-1])
df.low.Zhao2015I$Proportion <- factor(df.low.Zhao2015I$Proportion, levels = c("1%", "3%", "5%", "10%", "20%"))
df.low.Zhao2015I$Group <- factor(df.low.Zhao2015I$Group, levels = c("AIC", "BIC", "adjR2", "Cp"))
df.low.Zhao2015I$Perturbation <- factor(df.low.Zhao2015I$Perturbation, levels = c("10", "30", "50"))

# low-dimensional
# perturbation model IV (Masoud's model) 
df.low.Zhang2022IV <- as.data.frame(read.csv("Example I_Simulation I_dat_Zhang2022IV.csv")[,-c(1,5)])
df.low.Zhang2022IV$Proportion <- factor(df.low.Zhang2022IV$Proportion, levels = c("1%", "3%", "5%", "10%", "20%"))
df.low.Zhang2022IV$Group <- factor(df.low.Zhang2022IV$Group, levels = c("AIC", "BIC", "adjR2", "Cp"))
df.low.Zhang2022IV$Perturbation <- "IV"

# high-dimensional
# perturbation model I 
df.high.Zhao2015I <- as.data.frame(read.csv("Example I_Simulation II_all_dat_Zhao2015I.csv")[,-1])
df.high.Zhao2015I$Proportion <- factor(df.high.Zhao2015I$Proportion, levels = c("1%", "3%", "5%", "10%", "20%"))
df.high.Zhao2015I$Group <- factor(df.high.Zhao2015I$Group, levels = c("LASSO","SLASSO","ENET", "SCAD", "MCP"))
df.high.Zhao2015I$Perturbation <- factor(df.high.Zhao2015I$Perturbation, levels = c("10", "30", "50"))

# high-dimensional
# perturbation model IV (Masoud's model) 
df.high.Zhang2022IV <- as.data.frame(read.csv("Example I_Simulation II_dat_Zhang2022IV.csv")[,-c(1,5)])
df.high.Zhang2022IV$Proportion <- factor(df.high.Zhang2022IV$Proportion, levels = c("1%", "3%", "5%", "10%", "20%"))
df.high.Zhang2022IV$Group <- factor(df.high.Zhang2022IV$Group, levels = c("LASSO","SLASSO","ENET", "SCAD", "MCP"))
df.high.Zhang2022IV$Perturbation <- "IV"

#
dat.low <- rbind(df.low.Zhao2015I, df.low.Zhang2022IV)
dat.high <- rbind(df.high.Zhao2015I, df.high.Zhang2022IV)

count <- as.data.frame(read.csv("Example I_Simulation III_count.csv")[-c(1,10:51),-1])
count$Index <- c("1","2","3","4","5","6","...","100")
count$Index <- factor(count$Index, levels = c("1","2","3","4","5","6","...","100"))

################################################################################

jpeg(filename=paste0("Motivating_Example.jpeg"), width=2400, height=1100) 

font.size <- 45

txt10 <- expression(paste("I: ", kappa[paste(linear, ", ",bold(Y))], "=10")) 
txt30 <- expression(paste("I: ", kappa[paste(linear, ", ",bold(Y))], "=30")) 
txt50 <- expression(paste("I: ", kappa[paste(linear, ", ",bold(Y))], "=50"))
txtIV <- "IV"
txt.xlab <- expression(paste(zeta,":=(", n[infl],"/n)"%*%"100%")) 

fig1 <- ggplot(dat.low, aes(x=Proportion, y=Prob, group=interaction(Group, Perturbation))) +
  #geom_line(aes(linetype=1)) +
  geom_line(linetype = "dashed") + 
  geom_point(aes(shape=Perturbation, color=Group), size = 8)  +
  xlab(txt.xlab) + ylab("Empirical Probability") + 
  theme(axis.title.x = element_text(color="black", size=font.size, face="bold", vjust = -2),
        axis.title.y = element_text(color="black", size=font.size, face="bold", vjust = +3),
        axis.text.x = element_text(color="black",  size=font.size),
        axis.text.y = element_text(color="black",  size=font.size),
        #legend.title = element_blank(),  
        legend.text = element_text(size=font.size),
        legend.title = element_text(size=font.size, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = unit(c(3,3,3,3), "lines")) + 
  scale_y_continuous(breaks = seq(0, 1.0, 0.1)) + 
  coord_cartesian(ylim = c(0, 1.0))  + 
  scale_colour_discrete(name="Model\nSelectors", labels=c("AIC", "BIC", expression("Adjusted R"^"2"), expression("Mallow's C"["p"]))) + 
  scale_shape_discrete(name="Perturbation\nModels", labels=c(txt10,txt30,txt50,txtIV)) + 
  theme(legend.text.align = 0)

fig2 <- ggplot(dat.high, aes(x=Proportion, y=Prob, group=interaction(Group, Perturbation))) +
  geom_line(linetype = "dashed") +
  geom_point(aes(shape=Perturbation, color=Group), size = 8)  +
  xlab(txt.xlab) + ylab("Empirical Probability") + 
  theme(axis.title.x = element_text(color="black", size=font.size, face="bold", vjust = -2),
        axis.title.y = element_text(color="black", size=font.size, face="bold", vjust = +3),
        axis.text.x = element_text(color="black",  size=font.size),
        axis.text.y = element_text(color="black",  size=font.size),
        #legend.title = element_blank(),  
        legend.text = element_text(size=font.size),
        legend.title = element_text(size=font.size, face = "bold"),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        plot.margin = unit(c(3,3,3,3), "lines")) + 
  scale_y_continuous(breaks = seq(0, 1.0, 0.1)) + 
  coord_cartesian(ylim = c(0, 1.0)) + 
  scale_colour_discrete(name="Model\nSelectors", labels=c("LASSO","SLASSO","ENET", "SCAD", "MCP")) + 
  scale_shape_discrete(name="Perturbation\nModels", labels=c(txt10,txt30,txt50,txtIV)) + 
  theme(legend.text.align = 0)

final <- ggarrange(fig1, fig2,#fig3,  
                   labels = c("(a): low-dimensional", 
                              "(b): high-dimensional"),
                   ncol = 2, nrow = 1, vjust=-0.05, hjust=-0.60, common.legend = FALSE, font.label=list(color="black", size=42, face="bold"))

annotate_figure(final, top = text_grob("", color = "black", face = "bold", size = font.size))

dev.off()

#Empirical Probability of Selecting an Incorrect Submodel

#fig3 <- ggplot(count, aes(x=Index, y=Prob)) +
#  #geom_line(linetype="dashed") +
#  geom_point(size = 4, shape=23,fill="blue")  +
#  xlab("Indices of Observations") + ylab("") +
#  theme(axis.title.x = element_text(color="black", size=34, face="bold"),
#        axis.title.y = element_text(color="black", size=34, face="bold"),
#        #axis.text.x = element_text(color="black", size=25, angle = 45, hjust=1),
#        axis.text.x = element_text(color="black", size=34),
#        axis.text.y = element_text(color="black",  size=34),
#        #plot.title = element_text(color="black", size=27, hjust = 0.5, face="bold"),
#        plot.title = element_text(color="black", size=34, face="bold"),
#        panel.border = element_rect(color = "black", fill = NA, size = 1),
#        plot.margin = unit(c(3,3,3,3), "lines")) +
#  scale_y_continuous(breaks = seq(0, 1.0, 0.1)) + 
#  coord_cartesian(ylim = c(0, 1.0)) 
  
#  ,
#  legend.position = "bottom",
#  legend.box="vertical"  

#,
#legend.position = "bottom",
#legend.box="vertical"
