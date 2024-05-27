rm(list = ls())

##------------------------------------------------------------------------------
## Load libraries and source files
##------------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)
require(gridExtra)
library(ggpubr)
library(foreign)
library(xgboost)
library(fastglm)
library(randomForest)
library(grf)
library(doMC)
library(glmnet)
library(doParallel)
library(RColorBrewer)

# Call source files
source( "functions.R" ) 

##------------------------------------------------------------------------------
## Setup
##------------------------------------------------------------------------------
# Set seed
set.seed(123)

# Register cores for paralellization for lasso
num_cores = detectCores()
registerDoMC(cores = num_cores)
# Parallel backend
registerDoParallel(cores = num_cores)

# Color palette for figures
cb_colors = brewer.pal(n = 8, name = "Dark2") # discrete colorblind palette

##------------------------------------------------------------------------------
## True distribution
##------------------------------------------------------------------------------
# Approximate true distribution
df.true = dgp(1000000, 0.5)
ymin = min(df.true$y)
ymax = max(df.true$y)
# Locations for DTE, quantiles 0.1, ..., 0.9
vec.loc  = quantile(df.true$y, seq(0.1, 0.9, by=0.1))    
n.loc = length(vec.loc)
## Distributional treatment effect  (DTE)
dte.true = matrix(NA, n.loc, 1)
for(j in 1: n.loc ){
    v.loc   = vec.loc[j] ## location 
    dte.true[j]  = mean( (df.true$y[which(df.true$d==1)] < v.loc) ) -
              mean( (df.true$y[which(df.true$d==0)] < v.loc) ) 
      
}
# Save true DTE
saveRDS(dte.true, file=paste0("../result/dte.true.rds"))

##------------------------------------------------------------------------------
## Simulation setup
##------------------------------------------------------------------------------
# Number of simulations
n.sim  = 1000 
# Number of folds for cross-fitting
F = 5      
# Sample sizes
vec.n    = c(500, 1000, 5000)
# Treatment assignment probability
rho  = 0.5
# Gap length for PTE
h.pte   = 1  

##------------------------------------------------------------------------------
## Run simulation for DTE 
##------------------------------------------------------------------------------
start_time = Sys.time()
for(j in 1:length(vec.n)) { ## loop: sample size 
  n = vec.n[j]  ## sample size 
  # F-fold cross-fitting setup
  folds = sample(c(1:F), n, replace=TRUE)
  # Get results for each simulation run
  sim.output = foreach(i=1:n.sim) %dopar% run_simulation_dte(i)
}
end_time = Sys.time()
print(paste("Time spent:", end_time-start_time))

##------------------------------------------------------------------------------
## Run simulation for QTE 
##------------------------------------------------------------------------------
start_time = Sys.time()
for(j in 1:length(vec.n)) { ## loop: sample size 
  n = vec.n[j]  ## sample size 
  # F-fold cross-fitting setup
  folds = sample(c(1:F), n, replace=TRUE)
  # Get results for each simulation run
  sim.output = foreach(i=1:n.sim) %dopar% run_simulation_qte(i)
}
end_time = Sys.time()
print(paste("Time spent:", end_time-start_time))

##------------------------------------------------------------------------------
## Run simulation for sequence of DGPs
##------------------------------------------------------------------------------
start_time = Sys.time()
for (dgp_seq in 1:10){
  dgp_number = dgp_seq
  df.true = dgp(1000000, 0.5, dgp_number)
  ## True DTE
  vec.loc  = quantile(df.true$y, seq(0.1, 0.9, by=0.1))    ## locations for DTE
  n.loc = length(vec.loc)
  dte.true = matrix(NA, n.loc, 1)
  for(j in 1: n.loc ){
    v.loc   = vec.loc[j] ## location 
    dte.true[j]  = mean( (df.true$y[which(df.true$d==1)] < v.loc) ) -
      mean( (df.true$y[which(df.true$d==0)] < v.loc) ) 
    
  }
  saveRDS(dte.true, file=paste0("../result/", "dgp_seq", dgp_number, "_dte.true", ".rds"))
  
  n = 1000  ## sample size 
  # F-fold cross-fitting setup
  folds = sample(c(1:F), n, replace=TRUE)
  # Get results for each simulation run
  sim.output = foreach(i=1:n.sim) %dopar% run_simulation_sequence(i)
}
end_time = Sys.time()
print(paste("Time spent:", end_time-start_time))

