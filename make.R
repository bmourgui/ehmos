############################################%
# Run all the workflow
############################################%

#### Charge packages and functions ----
library(ggplot2)
library(magrittr)
####

# Load all functions
source(here::here("R", "analysis_functions.R"))
source(here::here("R", "graph_functions.R"))

#### Create folders to save results ----
dir.create("results")
dir.create(here::here("results", "figs"))
dir.create(here::here("results", "tables"))

#### File states ----
# Indicate if long-time-to-run files have changed
# Code below will be executed in consequence
simu01 <- "no"
simu02.glmm <- "ok"
simu02.ehmos <- "ok"
simu03 <- "no"
simu04 <- "no"

appli02 <- "no"

#### 1. Execute the workflow for simulation analysis ----
#### 1.1 - simulate the data ----
if (simu01 == "ok"){
  load(here::here("results", "simulated_data.RData"))
}else{
  source(here::here("analysis", "simu01_simulate-data.R"))
}


#### 1.2 - fit models
## GLMM
if (simu02.glmm == "ok"){
  load(here::here("results", "simu_out-glmm.RData"))
}else{
  source(here::here("analysis", "simu02_run-GLMM.R"))
}
## EHMOS
if (simu02.ehmos == "ok"){
  load(here::here("results", "simu_out-ehmos.RData"))
}else{
  source(here::here("analysis", "simu02_run-EHMOS.R"))
}

## 1.3. Summarise outputs
if (simu03 == "ok"){
  load(here::here("results", "simu_res-summary.RData"))
}else{
  load(here::here("results", "simu_out-ehmos.RData"))
  load(here::here("results", "simu_out-glmm.RData"))
  source(here::here("analysis", "simu03_summarise-model-outputs.R"))
}

#### 1.4 - Compute performance metrics ----
if (simu04 == "ok"){
  load(here::here("results", "simu_res-perf.RData"))
}else{
  source(here::here("analysis", "simu04_compute-performance-metrics.R"))
}

source(here::here("analysis", "simu05_make-figs&tabs-results.R"))


#### 2. Execute workflow for the case study ----
#### .. 2.1. Data formatting ----
source(here::here("analysis", "appli01_formatting-data.R"))

#### .. 2.2. fit models ----
if (appli02 == "ok"){
  load(here::here("results", "appli_out-glmm.RData"))
  load(here::here("results", "appli_out-ehmos.RData"))
}else{
  source(here::here("analysis", "appli02_run-models.R"))
}

#### .. 2.3. Summarise model outputs ----
source(here::here("analysis" , "appli03_summarise-model-outputs.R"))

#### .. 2.4. Results analysis ----
source(here::here("analysis" , "appli04_results.R"))


#knitr::knit(input = here::here("redaction", "manuscript.Rmd"))
#rmarkdown::render(input = here::here("redaction", "manuscript.Rmd"), rmarkdown::pdf_document())
