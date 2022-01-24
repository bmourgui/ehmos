############################################%
# Calculate performance metrics
#
# last modification: 09/12/21 (edit annotations)
# bastien.mourguiart@gmail.com
#
# simu04_compute-performance-metrics.R
# depends on: simu03_summarise-model-outputs.R
#
# Script to calculate performance metrics based
# on summarized model outputs
############################################%

## Change scenario names, to scope with final ones
full.res.shift$scenario <- as.factor(full.res.shift$scenario)
levels(full.res.shift$scenario) <- c("A2xB1xC1","A2xB1xC2","A2xB2xC1","A2xB2xC2","A1xB1xC1","A1xB1xC2","A1xB2xC1","A1xB2xC2")

## Add contrain.true column that will be used to calculate
## interval coverage
full.res.shift %>%
  dplyr::mutate("contain.true" = contains(lwr.shift,
                                          upr.shift,
                                          true.shift)) -> full.res.shift

## Add interval score
full.res.shift$IS <- NA
for (i in 1:nrow(full.res.shift)){
  full.res.shift[i,]$IS <- ISfunction(lwr = full.res.shift[i,]$lwr.shift,
                                      upr = full.res.shift[i,]$upr.shift,
                                      true = full.res.shift[i,]$true.shift)
}


## Performance at the species level
full.res.shift %>%
  dplyr::group_by(model,
                  scenario,
                  sp,
                  T_opti,
                  T_shape
                  ) %>%
  dplyr::summarise("RMSE"=sqrt(sum((estim.shift-true.shift)^2)/R),
                   "SE"=sqrt(sum((estim.shift-mean(estim.shift))^2)/(R-1)),
                   "bias"=mean(abs(estim.shift)-abs(true.shift)), # take the absolute value to see bias on the magnitude of the shift -> a negative bias will involve an underestimation of shift
                   "coverage"=mean(sum(contain.true)/dplyr::n()),
                   "IS" = mean(IS)
                   ) %>%
  dplyr::mutate("sampling" = substr(scenario, 1, 2),
                "marginality" = substr(scenario, 4, 5),
                "specialization" = substr(scenario, 7, 8)
                ) -> perf_sp

## Overall performance: average at the model level
perf_sp %>%
  average_metrics(model) -> overall_perf

## Community-level performance: average at the scenario level
perf_sp %>%
  average_metrics(model, scenario) -> perf_scenario

## Optimum-type level performance: average species according to their marginality type
perf_sp %>%
  average_metrics(model, T_opti) -> perf_Topti

save(perf_sp,
     overall_perf,
     perf_scenario,
     perf_Topti,
     file = here::here("results", "simu_res-perf.RData"))
