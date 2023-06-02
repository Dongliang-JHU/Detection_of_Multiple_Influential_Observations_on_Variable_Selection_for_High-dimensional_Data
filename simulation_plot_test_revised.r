
library(ggplot2)
library(stringr)
library(gridExtra) 
library(lemon) #(grid_arrange_shared_legend)  
library(dplyr)
library(ggpubr)

################################################################################

development.dir <- "/Users/dongliangzhang/Library/CloudStorage/Dropbox/Research/Detection of Influential Observations on Variable Selection/Dongliang Zhang_Influential Diagnostics/Code/Development/" 

source.dir <- "/Users/dongliangzhang/Library/CloudStorage/Dropbox/Research/Detection of Influential Observations on Variable Selection/Dongliang Zhang_Influential Diagnostics/Code/Simulation/"

output.dir <- "/Users/dongliangzhang/Library/CloudStorage/Dropbox/Research/Detection of Influential Observations on Variable Selection/Dongliang Zhang_Influential Diagnostics/Code/Output/"

output.dir <- development.dir

source(paste0(development.dir,"simulation_plot_revised.r"))

################################################################################
#Linear Model: Perturbation Model I 

sim.plot(source.dir, output.dir, X.type="IID", perturb.type="Zhao2015I", y.perturb=5, x.perturb=0)

sim.plot(source.dir, output.dir, X.type="IID", perturb.type="Zhao2015I", y.perturb=10, x.perturb=0)

sim.plot(source.dir, output.dir, X.type="AR5", perturb.type="Zhao2015I", y.perturb=5, x.perturb=0)

sim.plot(source.dir, output.dir, X.type="AR5", perturb.type="Zhao2015I", y.perturb=10, x.perturb=0)

sim.plot(source.dir, output.dir, X.type="AR8", perturb.type="Zhao2015I", y.perturb=5, x.perturb=0)

sim.plot(source.dir, output.dir, X.type="AR8", perturb.type="Zhao2015I", y.perturb=10, x.perturb=0)

################################################################################

#Linear Model: Perturbation Model II  

sim.plot(source.dir, output.dir, X.type="IID", perturb.type="Zhang2022II", y.perturb=0, x.perturb=5)

sim.plot(source.dir, output.dir, X.type="IID", perturb.type="Zhang2022II", y.perturb=0, x.perturb=10)

sim.plot(source.dir, output.dir, X.type="AR5", perturb.type="Zhang2022II", y.perturb=0, x.perturb=5)

sim.plot(source.dir, output.dir, X.type="AR5", perturb.type="Zhang2022II", y.perturb=0, x.perturb=10)

sim.plot(source.dir, output.dir, X.type="AR8", perturb.type="Zhang2022II", y.perturb=0, x.perturb=5)

sim.plot(source.dir, output.dir, X.type="AR8", perturb.type="Zhang2022II", y.perturb=0, x.perturb=10)

################################################################################

#Linear Model: Perturbation Model III  

sim.plot(source.dir, output.dir, X.type="IID", perturb.type="Zhao2015III", y.perturb=5, x.perturb=5)

sim.plot(source.dir, output.dir, X.type="IID", perturb.type="Zhao2015III", y.perturb=10, x.perturb=10)

sim.plot(source.dir, output.dir, X.type="AR5", perturb.type="Zhao2015III", y.perturb=5, x.perturb=5)

sim.plot(source.dir, output.dir, X.type="AR5", perturb.type="Zhao2015III", y.perturb=10, x.perturb=10)

sim.plot(source.dir, output.dir, X.type="AR8", perturb.type="Zhao2015III", y.perturb=5, x.perturb=5)

sim.plot(source.dir, output.dir, X.type="AR8", perturb.type="Zhao2015III", y.perturb=10, x.perturb=10)

################################################################################

#Poisson Model: Perturbation Model I  

sim.plot(source.dir, output.dir, X.type="IID", perturb.type="PoissonI", y.perturb=1, x.perturb=0)

sim.plot(source.dir, output.dir, X.type="AR5", perturb.type="PoissonI", y.perturb=1, x.perturb=0)

sim.plot(source.dir, output.dir, X.type="AR8", perturb.type="PoissonI", y.perturb=1, x.perturb=0)

################################################################################

#Poisson Model: Perturbation Model II  

sim.plot(source.dir, output.dir, X.type="IID", perturb.type="Poisson2023II", y.perturb=0, x.perturb=10)

################################################################################

#Poisson Model: Perturbation Model III  

sim.plot(source.dir, output.dir, X.type="IID", perturb.type="Poisson2023III", y.perturb=10, x.perturb=10)



