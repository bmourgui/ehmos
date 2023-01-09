# Explicit Modelling of Optimum Shifts

This repository contains the data and code to reproduce analysis presented in the paper:

> B. Mourguiart, B. Liquet, K. Mengersen, T. Couturier, J. Mansons, Y. Braud, A. Besnard. A new method to explicitly estimate the shift of optimum along gradients in multispecies studies. *Journal of Biogeography*, (in press).

The paper introduces a new formulation of a Bayesian hierarchical linear model that explicitly estimates optimum shifts for multiple species having symmetrical response curves. This new formulation, called Explicit Hierarchical Model of Optimum Shifts (EHMOS), is compared to a mean comparison method and a Bayesian generalized linear mixed model (GLMM) using simulated and real datasets.

# Home directory

Important contents:
- This README
- ehoms_github.Rproj: the R project file. Download the entire repository and open this file in R studio to run all R scripts with the correct relative paths.
- make.R: the R script workflow. Run this script to run all the scripts in the correct order to reproduce all analysis presented in the paper. 

# R
This folder contains scripts coding the home-made functions required for analysis.

# data
This folder contains occurrence data used in the application study on Orthoptera species. One file for each survey.

# analysis
This folder contains all R scripts needed to reproduce the results of the simulation study (scripts begining with *simu*), the application study (scripts begining with *appli*) and the supplementary analysis (scripts begining with *supp*). The numbers indicated the order in which scripts should be run.
- simu
  - simu01_simulate-data.R: simulate the occurrence data for the 20 virtual species studied in the simulation study.
  - simu02-run-EHMOS.R: code to fit the Explicit Hierarchical Model of Optimum Shifts (EHMOS) to the simulated data. Instead of running this script which can take several days to run, the model's outputs (MCMC samples) used to produce the results in the paper can be download at: . 
  - simu02-run-GLMM.R: code to fit the GLMM to the simulated data. Instead of running this script which can take several days to run, the model's outputs (MCMC samples) used to produce the results in the paper can be download at: . 
  - simu03_summarise-model-outputs.R: transforms models' outputs (MCMC samples) in a summary dataframe containing optimum shift estimates for each species, models, scenarios and replications. 
  - simu04_compute-performance-metrics.R: computes performance metrics (e.g., bias)  for each species, models, scenarios and replications. 
  - simu05_make-figs&tabs-results.R: creates figures and tables presented in the main paper and supplementary materials. 
