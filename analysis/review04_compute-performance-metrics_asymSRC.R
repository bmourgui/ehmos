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
                  curve.type1,
                  curve.type2
                  ) %>%
  dplyr::summarise("RMSE"=sqrt(sum((estim.shift-true.shift)^2)/R),
                   "SE"=sqrt(sum((estim.shift-mean(estim.shift))^2)/(R-1)),
                   "bias"=mean(abs(estim.shift)-abs(true.shift)), # take the absolute value to see bias on the magnitude of the shift -> a negative bias will involve an underestimation of shift
                   "coverage"=mean(sum(contain.true)/dplyr::n()),
                   "IS" = mean(IS)
                   ) %>%
  dplyr::mutate("sptype" = paste(curve.type1, curve.type2)) -> perf_sp

## Overall performance: average at the model level
perf_sp %>%
  average_metrics(model) -> overall_perf

## Community-level performance: average at the scenario level
perf_sp %>%
  average_metrics(model, scenario) -> perf_scenario

## Optimum-type level performance: average species according to their marginality type
perf_sp %>%
  average_metrics(model, sptype) -> perf_sptype

save(perf_sp,
     overall_perf,
     perf_scenario,
     perf_Topti,
     file = here::here("results", "simu_asym_res-perf.RData"))
