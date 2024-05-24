
## Estimating Distributional Treatment Effects in Randomized Experiments: Machine Learning for Variance Reduction

### Folders 

1. `data` folder includes the experiment data from Ferraro & Price
(2013)

2. `result` folder is empty and all intermediate csv files and figures
will be saved in the folder

### Files 

1. `functions.R` file includes all necessary functions

2. `run_simulation.R` includes code to run the Monte Carlo simulations and save results as rds files

3. `compute_stats.R` includes code to calculate measures (e.g. bias, RMSE) from the saved simulation results (rds files) and save them as csv files

4. `plot_figures.R` includes code to load the csv files and plot figures for the simulation study

3. `experiment_water_consumption.R` includes code to replicate the analysis of experimental data from Ferraro & Price (2013)

### Instructions

1.  Install all necessary packages in R
2.  To replicate the results from the
    Monte Carlo simulation, run the files in the following order: `run_simulation.R` -> `compute_stats.R` -> `plot_figures.R`
The outputs will be figures appeared in Figures 1, 3 and 4 in the paper. 
    in the paper
3.  Run `experiment_water_consumption.R`file to replicate the results from the experimental analysis. The output will be figures appeared in
    Figure 2 in the paper
