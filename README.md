
## Estimating Distributional Treatment Effects in Randomized Experiments: Machine Learning for Variance Reduction
This repository contains code to replicate the experimental results from "Estimating Distributional Treatment Effects in Randomized Experiments: Machine Learning for Variance Reduction."

### Folders 

1. `data` folder includes files to create dataset used for empirical application from [Ferraro & Price (2013)](https://direct.mit.edu/rest/article-abstract/95/1/64/58053/Using-Nonpecuniary-Strategies-to-Influence)

2. `experiment` folder contains all R files used for analysis

### Experiment Files 

1. `functions.R` file includes all necessary functions

2. `run_simulation.R` includes code to run the Monte Carlo simulations and saves results as .rds files

3. `compute_stats.R` includes code to calculate evaluation metrics (e.g. bias, RMSE) from the saved simulation results (.rds files) and saves them as .csv files

4. `plot_figures.R` includes code to load the .csv files and plot figures for the simulation study

3. `experiment_water_consumption.R` includes code to replicate the analysis of experimental data from Ferraro & Price (2013)

### Instructions

1.  Install all necessary packages in R
2.  To replicate the results from the
    Monte Carlo simulation, run the files in the following order: (1) `run_simulation.R`, (2) `compute_stats.R`, (3) `plot_figures.R`. The outputs will be figures appeared in Figures 1, 3 and 4 in the paper. 
3.  Run `experiment_water_consumption.R` to replicate the results from the water consumption experiment. The output will be figures appeared in
    Figure 2 in the paper.
    
### R version and attached packages
- R version 4.3.1

- `RColorBrewer_1.1-3` `ggpubr_0.6.0`       `fastglm_0.0.3`      `bigmemory_4.6.1`    `xgboost_1.7.5.1`    `foreign_0.8-84`     `ggplot2_3.4.3` `dplyr_1.1.2`  `doParallel_1.0.17`    `glmnet_4.1-8`         `Matrix_1.6-1.1`       `doMC_1.3.8`           `iterators_1.0.14`     `foreach_1.5.2`        `grf_2.3.1` `randomForest_4.7-1.1` `gridExtra_2.3`        `tidyr_1.3.0` 
