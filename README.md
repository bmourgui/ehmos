# Explicit Modelling of Optimum Shifts

This repository contains the data and code to reproduce analysis presented in the paper:

> B. Mourguiart, B. Liquet, K. Mengersen, T. Couturier, J. Mansons, Y. Braud, A. Besnard. A new method to explicitly estimate the shift of optimum along gradients in multispecies studies. *Journal of Biogeography*, (in press), doi: 10.1111/jbi.14570.

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
This folder contains all R scripts needed to reproduce the results of the simulation study (scripts begining with *simu*), the application study (scripts begining with *appli*) and the supplementary analysis about assumption of symmetrical responses (scripts begining with *supp*). The numbers indicated the order in which scripts should be run.

- simu
  - simu01_simulate-data.R: simulate the occurrence data for the 20 virtual species studied in the simulation study.
  - simu02-run-EHMOS.R: code to fit the Explicit Hierarchical Model of Optimum Shifts (EHMOS) to the simulated data. Instead of running this script which can take several days to run, the model's outputs (MCMC samples saved in *simu_out-ehmos.RData*) used to produce the results in the paper can be download at: https://doi.org/10.5061/dryad.0gb5mkm59 (file simu_out-ehmos.RData). 
  - simu02-run-GLMM.R: code to fit the GLMM to the simulated data. Instead of running this script which can take several days to run, the model's outputs (MCMC samples saved in *simu_out-glmm.RData*) used to produce the results in the paper can be download at: https://doi.org/10.5061/dryad.0gb5mkm59 (file simu_out-glmm.RData). 
  - simu03_summarise-model-outputs.R: transforms models' outputs (MCMC samples) in a summary dataframe containing optimum shift estimates for each species, models, scenarios and replications. 
  - simu04_compute-performance-metrics.R: computes performance metrics (e.g., bias)  for each species, models, scenarios and replications. 
  - simu05_make-figs&tabs-results.R: creates figures and tables presented in the main paper and supplementary materials. 

- appli
  - appli01_formatting-data.R: transform raw datasets store in the data folder into a matrix used to fit the models.
  - appli02_run-models.R: fit the two Bayesian models and the mean comparison method to the data created in appli01_formatting-data.R.
  - appli03_summarise-model-outputs.R: create a data frame containing summary statistics of estimates (e.g., mean shift estimates) obtained with the three methods.
  - appli04_results.R: create main results of the application study.

- supp
  - supp01_simulate-data.R: simulate the occurrence data for 9 virtual species with both symmetrical and asymmetrical response curves.
  - supp02-run-EHMOS.R: code to fit the Explicit Hierarchical Model of Optimum Shifts (EHMOS) to the simulated data. 
  - supp02-run-GLMM.R: code to fit the GLMM to the simulated data. 
  - supp03_summarise-model-outputs.R: transforms models' outputs (MCMC samples) in a summary dataframe containing optimum shift estimates for each species, models and replications. 
  - supp04_compute-performance-metrics.R: computes performance metrics (e.g., bias)  for each species, models and replications. 
  - supp05_make-figs&tabs-results.R: creates figures and tables presented in supplementary materials 4. 
