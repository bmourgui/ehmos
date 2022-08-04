############################################%
# Summarise model outputs
#
# last modification: 09/12/21 (edit annotations)
# bastien.mourguiart@gmail.com
#
# simu03_summarise-model-outputs.R
#
# Script to summarise model outputs
# Format data for visualization
############################################%

library(magrittr)

source(here::here("R/analysis_functions.R"))

load(here::here("results/simulated_data_asymSRC.RData"))
load(here::here("results/simu_asym_out-glmm.RData"))
load(here::here("results/simu_asym_out-ehmos.RData"))


data.sp$sp <- paste0("sp", 1:9)

# Set indices
N <-  dim(Z)[2] # number of speceis
J <- dim(Z)[1] # number of sites
S <- dim(Z)[3] # number of scenarios
R <- 30  #nb of repli
M <- 3 #nb of methods

#### 1. Summarise estimate outputs ####
# Create data frame to fill
full.res.shift <- data.frame("model" = rep(c("ehmos", "t.test", "glmm"),
                                           each = N
                                           ),
                             "scenario" = rep(c("A1", "A2"),
                                              each = N*R*M
                                              ),
                             "replication" = rep(1:R, each = N*M),
                             "sp" = paste0("sp", 1:N),
                             "converge" = NA, # boolean vector telling if all parameters converged
                             "time" = NA, # computation time in seconds
                             "estim.shift" = NA, # shift estimates
                             "lwr.shift" = NA, # lower boundary of 95% credible interval
                             "upr.shift" = NA # upper boundary of 95% credible interval)
                             ) %>%
  dplyr::inner_join(data.sp,
    by = c("sp")
  ) %>%
  dplyr::rename(true.shift = shift)  # do not distinguish 'edge down' (species with optima close to the bottom sampling edge) from 'edge up'



for (s in 1:S) {
  for (r in 1:R) {

    name.ehmos <- paste0("out.sim_asym_A", s, "_EHMOS", r)
    name.glmm <- paste0("out.sim_asym_A", s, "_glmm", r)
    out <- res.ehmos_asym[[name.ehmos]]
    out.glmm <- res.glmm_asym[[name.glmm]]

    x <- get(paste0("X_A",s)) # extract appropriate elevation vector
    z <- Z[ , , , s, r]

    # Calculate first optima
    opti1.ehmos <- unscale(param = out$sims.list$opt, # directly estimated
                           x = x)
    opti1.glmm <- unscale(param = -out.glmm$sims.list$b1/(2*out.glmm$sims.list$b2), # derive optima from model parameters
                          x = x)

    # Calculate second optima
    opti2.glmm <- unscale(param = -(out.glmm$sims.list$b1 + out.glmm$sims.list$b4) / (2*(out.glmm$sims.list$b2 + out.glmm$sims.list$b5)),
                          x = x)

    # Derive shifts
    shift.ehmos <- out$sims.list$shift*sd(x)
    shift.glmm <- opti2.glmm - opti1.glmm

    # Fill the data frame
    full.res.shift[full.res.shift$scenario == paste0("A", s) & full.res.shift$replication==r & full.res.shift$model == "ehmos", "converge"] <- sum(lapply(out$Rhat, function(x)sum(x>1.1))>0)
    full.res.shift[full.res.shift$scenario == paste0("A", s) & full.res.shift$replication==r & full.res.shift$model == "t.test", "converge"] <- 0
    full.res.shift[full.res.shift$scenario == paste0("A", s) & full.res.shift$replication==r & full.res.shift$model == "glmm", "converge"] <- sum(lapply(out.glmm$Rhat, function(x)sum(x>1.1))>0)

    full.res.shift[full.res.shift$scenario == paste0("A", s) & full.res.shift$replication==r & full.res.shift$model == "ehmos", "time"] <- out$mcmc.info$elapsed.mins
    full.res.shift[full.res.shift$scenario == paste0("A", s) & full.res.shift$replication==r & full.res.shift$model == "t.test", "time"] <- 0
    full.res.shift[full.res.shift$scenario == paste0("A", s) & full.res.shift$replication==r & full.res.shift$model == "glmm", "time"] <- out.glmm$mcmc.info$elapsed.mins

    full.res.shift[full.res.shift$scenario == paste0("A", s) & full.res.shift$replication==r & full.res.shift$model == "ehmos", c("estim.shift", "lwr.shift", "upr.shift")] <- t(apply(shift.ehmos, 2, function(y)quantile(y, c(0.5,0.025,0.975))))
    full.res.shift[full.res.shift$scenario == paste0("A", s) & full.res.shift$replication==r & full.res.shift$model == "t.test", c("lwr.shift", "estim.shift", "upr.shift")] <- t(apply(z, 2, function(y)mean.method2(z=y, env=x)))[,1:3]
    full.res.shift[full.res.shift$scenario == paste0("A", s) & full.res.shift$replication==r & full.res.shift$model == "glmm", c("estim.shift", "lwr.shift", "upr.shift")] <- t(apply(shift.glmm, 2, function(y)quantile(y, c(0.5,0.025,0.975))))
  }
}

## Add the number of presence per species x scenario x replication
apply(Z, MARGIN = c(2,3,4, 5), FUN = sum) %>%
  reshape2::melt() %>%
  dplyr::rename(sp = Var1,
                period = Var2,
                scenario = Var3,
                replication = Var4,
                n = value) %>%
  dplyr::mutate("replication" = as.numeric(substr(replication, 6, length(replication)))
  ) %>%
  tidyr::spread(period, n) %>%
  dplyr::rename(n1 = p1,
                n2 = p2) %>%
  dplyr::inner_join(full.res.shift,
                    by = c("scenario",
                           "replication",
                           "sp")) -> full.res.shift

save(full.res.shift, file = here::here("results", "simu_asym_res-summary.RData"))

