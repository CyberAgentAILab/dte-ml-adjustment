rm(list = ls())

##------------------------------------------------------------------------------
## Load libraries and source files
##------------------------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(foreign)
library(xgboost) 
library(fastglm)
library(randomForest)
library(grf)
library(doMC) 
library(glmnet)
library(ggpubr) 
library(doParallel)

##------------------------------------------------------------------------------
## Setup
##------------------------------------------------------------------------------
vec.n = c(500, 1000, 5000)
vec.loc = seq(0.1, 0.9, 0.1)

##------------------------------------------------------------------------------
## DTE simulation
##------------------------------------------------------------------------------
for(j in 1:length(vec.n)) { ## loop: sample size 
  
  n = vec.n[j]  ## sample size 
  
  # Load true DTE
  dte.true = readRDS(paste0("./result/dte.true.rds"))
  
  # Set working directory
  setwd(paste0("./result/dte_n", n, "/"))
  # List all RDS files in the directory
  rds_files = list.files(pattern = "\\.rds$")
  # Initialize an empty list to store loaded data frames
  results = list()
  # Load each RDS file and store it in the list
  for (file in rds_files) {
    results_eachround = readRDS(file)
    results[[file]] = results_eachround
  }
  
  # Number of simulations
  n.sim = length(rds_files)
  print(n.sim)

  # RMSE
  rmse0 = sqrt(Reduce('+', lapply(results, function(mat)(mat$res.ols$dte$dte-dte.true)^2))/n.sim)
  rmse1 = sqrt(Reduce('+', lapply(results, function(mat)(mat$res.ols$dte.ra$dte-dte.true)^2))/n.sim)
  rmse2 = sqrt(Reduce('+', lapply(results, function(mat)(mat$res.lasso$dte.ra$dte-dte.true)^2))/n.sim)
  
  rmse  = cbind(n, vec.loc, rmse0, rmse1, rmse2)
  colnames(rmse) = c( "n","location","Simple","OLS", "Lasso")
  
  # Bias
  bias0 = Reduce('+', lapply(results, function(mat)(mat$res.ols$dte$dte-dte.true)))/n.sim
  bias1 = Reduce('+', lapply(results, function(mat)(mat$res.ols$dte.ra$dte-dte.true)))/n.sim
  bias2 = Reduce('+', lapply(results, function(mat)(mat$res.lasso$dte.ra$dte-dte.true)))/n.sim
  
  bias = cbind(n, vec.loc, bias0, bias1, bias2)
  colnames(bias) = c( "n","location","Simple","OLS", "Lasso")
  
  # RMSE ratio
  ratio = cbind(n, vec.loc, 100*(1-rmse1/rmse0), 100*(1-rmse2/rmse0))
  colnames(ratio) = c( "n","location","OLS", "Lasso")
  
  # Save results as csv files
  v.info = paste0("dte_n", n)
  write.csv(rmse, file=paste0(v.info, "_rmse.csv"), row.names=FALSE)
  write.csv(ratio, file=paste0(v.info, "_ratio.csv"), row.names=FALSE)
  write.csv(bias, file=paste0(v.info, "_bias.csv"), row.names=FALSE )
  write.csv(dte.true, file = paste0(v.info, "_dte.true.csv"), row.names=FALSE)
}

##------------------------------------------------------------------------------
## Simulation for QTE
##------------------------------------------------------------------------------
for(j in 1:length(vec.n)) { ## loop: sample size 
  
  n = vec.n[j]  ## sample size 
  
  # Set working directory
  setwd(paste0("./qte_n", n, "/"))
  # List all RDS files in the directory
  rds_files = list.files(pattern = "\\.rds$")
  # Initialize an empty list to store loaded data frames
  results = list()
  # Load each RDS file and store it in the list
  for (file in rds_files) {
    results_eachround = readRDS(file)
    results[[file]] = results_eachround
  }
  
  # Number of simulations
  n.sim = length(rds_files)
  
  # RMSE
  rmse0 = sqrt(Reduce('+', lapply(results, function(mat)(mat$res.ols$qte$qte-1)^2))/n.sim)
  rmse1 = sqrt(Reduce('+', lapply(results, function(mat)(mat$res.ols$qte.ra$qte-1)^2))/n.sim)
  rmse2 = sqrt(Reduce('+', lapply(results, function(mat)(mat$res.lasso$qte.ra$qte-1)^2))/n.sim)
  
  rmse  = cbind(n, vec.loc, rmse0, rmse1, rmse2)
  colnames(rmse) = c( "n","location","Simple","OLS", "Lasso")
  
  # Bias
  bias0 = Reduce('+', lapply(results, function(mat)(mat$res.ols$qte$qte-1)))/n.sim
  bias1 = Reduce('+', lapply(results, function(mat)(mat$res.ols$qte.ra$qte-1)))/n.sim
  bias2 = Reduce('+', lapply(results, function(mat)(mat$res.lasso$qte.ra$qte-1)))/n.sim
  
  bias  = cbind(n, vec.loc, bias0, bias1, bias2)
  colnames(bias) = c( "n","location","Simple","OLS", "Lasso")
  
  # RMSE ratio
  ratio = cbind(n, vec.loc, 100*(1-rmse1/rmse0), 100*(1-rmse2/rmse0))
  colnames(ratio) = c( "n","location","OLS", "Lasso")
  
  # Save results as csv files
  v.info = paste0("qte_n", n)
  write.csv(rmse, file=paste0( v.info, "_rmse.csv"), row.names=FALSE )
  write.csv(bias, file=paste0(v.info, "_bias.csv"), row.names=FALSE )
  write.csv(ratio, file=paste0(v.info, "_ratio.csv"), row.names=FALSE)
}

##------------------------------------------------------------------------------
## Simulation for sequence of DGPs
##------------------------------------------------------------------------------
for (dgp_seq in 1:10){
  
  n = 1000
  # Load true DTE
  dte.true = readRDS(paste0("./result/dgp_seq", dgp_number, "_dte.true.rds"))
  
  # Set working directory
  setwd(paste0("./result/dgp_seq", dgp_number, "/"))
  # List all RDS files in the directory
  rds_files = list.files(pattern = "\\.rds$")
  # Initialize an empty list to store loaded data frames
  results = list()
  # Load each RDS file and store it in the list
  for (file in rds_files) {
    results_eachround = readRDS(file)
    results[[file]] = results_eachround
  }
  
  # Number of simulations
  n.sim = length(rds_files)
  
  # RMSE
  rmse0 = sqrt(Reduce('+', lapply(results, function(mat)(mat$res.ols$dte$dte-dte.true)^2))/n.sim)
  rmse1 = sqrt(Reduce('+', lapply(results, function(mat)(mat$res.ols$dte.ra$dte-dte.true)^2))/n.sim)
  rmse2 = sqrt(Reduce('+', lapply(results, function(mat)(mat$res.lasso$dte.ra$dte-dte.true)^2))/n.sim)
  
  rmse  = round(cbind(n, vec.loc, rmse0, rmse1, rmse2), 2)
  colnames(rmse) = c( "n","location","Simple","OLS", "Lasso")
  
  # RMSE ratio
  ratio = round(cbind(n, vec.loc, 100*(1-rmse1/rmse0), 100*(1-rmse2/rmse0)),2)
  colnames(ratio) = c( "n","location","OLS", "Lasso")
  
  # Save results as csv files
  v.info = paste0("dgp_seq", dgp_number)
  write.csv(rmse, file=paste0( v.info, "_rmse.csv"), row.names=FALSE )
  write.csv(ratio, file=paste0( v.info, "_ratio.csv"), row.names=FALSE )
  write.csv(dte.true, file = paste0( v.info, "_dte.true.csv"), row.names=FALSE )
  
}




