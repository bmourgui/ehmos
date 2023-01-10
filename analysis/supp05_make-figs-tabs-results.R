
library(ggplot2)
library(magrittr)

load(here::here("results/simulated_data_asymSRC.RData"))
load(here::here("results/simu_asym_out-glmm.RData"))
load(here::here("results/simu_asym_out-ehmos.RData"))
load(here::here("results", "simu_asym_res-summary.RData"))

R <- 30
N <- 9
J <- 300
S <- 2

#### 1. Figure 1: Bias in species shift estimates ----
## Bias
full.res.shift %>%
  dplyr::mutate("bias" = (estim.shift - true.shift)) %>%
  dplyr::group_by(model, scenario, sp) %>%
  dplyr::summarise_at("bias", c("mean" = mean, "sd" = sd)) -> summary_bias

full.res.shift %>%
  dplyr::mutate("bias" = (estim.shift - true.shift)) %>%
  ggplot(aes(x = sp, y = bias)) +
  geom_point(color = alpha("grey20", 0.4)) +
  facet_grid(scenario ~ model,
             labeller = labeller(.rows = c("A1" = "A2", "A2" = "A1"))) + # correct past error of confounding names when creating data
  geom_point(data = summary_bias,
             aes(x = sp, y = mean), 
             color = "red",
             size = 2) +
  geom_segment(data = summary_bias,
               aes(x = sp, 
                   xend = sp, 
                   y = mean - sd,
                   yend = mean + sd), 
               color = "red") +
  ylab("Bias (m)") +
  theme_bw()

## Relative bias
full.res.shift %>%
  dplyr::mutate("bias" = (estim.shift - true.shift)/true.shift) %>%
  dplyr::group_by(model, scenario, sp) %>%
  dplyr::summarise_at("bias", c("mean" = mean, "sd" = sd)) -> summary_relative_bias

full.res.shift %>%
  dplyr::mutate("bias" = (estim.shift - true.shift)/true.shift) %>%
  ggplot(aes(x = sp, y = bias)) +
  geom_point(color = alpha("grey20", 0.4)) +
  facet_grid(scenario ~ model,
             labeller = labeller(.rows = c("A1" = "A2", "A2" = "A1"))) +
  geom_point(data = summary_relative_bias,
             aes(x = sp, y = mean), 
               color = "red",
             size = 2) +
  geom_segment(data = summary_relative_bias,
               aes(x = sp, 
                   xend = sp, 
                   y = mean - sd,
                   yend = mean + sd), 
               color = "red") +
  ylab("Relative bias") +
  theme_bw()

#### 1.2. Table of Bias results ####
summary_bias %>%
  dplyr::ungroup() %>%
  dplyr::inner_join(full.res.shift[, c("sp", "true.shift")], by = "sp") %>%
  dplyr::distinct() %>%
  dplyr::mutate("bias" = paste0(round(mean), " (", round(sd), ")")) %>%
  dplyr::select(-c('mean', "sd")) %>%
  tidyr::pivot_wider(names_from = model, values_from = bias) %>%
  dplyr::rename("Scenario" = "scenario", 
                "Species" = "sp", 
                "True shift" = "true.shift", 
                "EHMOS" = "ehmos", 
                "cGLMM" = "glmm", 
                "t-test" = "t.test")  %>%
  dplyr::select(-c("Scenario")) %>%
  kableExtra::kable(format = "latex",
                    caption = "Average of bias in optimum shift estimates for the two simulated scenarios after 30 replications. Numbers in brackets represent standard deviations. ",
                    booktabs = TRUE,
                    longtable = TRUE) %>%
  kableExtra::kable_styling(latex_options = "HOLD_position") %>%
  kableExtra::pack_rows("Sampling scenario A1", 1, 9) %>%
  kableExtra::pack_rows("Sampling scenario A2", 10, 18)


#### 2. Goodness-of-fit ####
out <- res.ehmos_asym[[1]]
p.fit <- out$sims.list$p.fit
p.fitnew <- out$sims.list$p.fitnew
bp <- sum(p.fit > p.fitnew)/length(p.fit)
plot(p.fit, p.fitnew,
     xlim = range(c(p.fit, p.fitnew)),
     ylim = range(c(p.fit, p.fitnew)))
abline(a = 0, b = 1)
mtext(paste0("Bayesian p-value = ", round(bp, 2)), adj = 0)

par(mfrow = c(3, 3))
for (i in 1:N){ # N, the number of species
  p.fit <- out$sims.list$p.fit.sp[,i]
  p.fitnew <- out$sims.list$p.fitnew.sp[,i]
  bp <- sum(p.fit > p.fitnew)/length(p.fit)
  plot(p.fit, p.fitnew,
       xlim = range(c(p.fit, p.fitnew)),
       ylim = range(c(p.fit, p.fitnew)),
       main = paste0("Species ", i))
  abline(a = 0, b = 1)
  mtext(paste0("Bayesian p-value = ", round(bp, 2)), adj = 0)
}

# count the number of replications where Bayes p-value indicates lack of fit
# for ehmos
bp_ehmos <- vector("numeric", length = R)
bp.sp_ehmos <- matrix(0, nrow = N, ncol = R)
for (r in 1:R){
  # general
  out_ehmos <- res.ehmos_asym[[r]]
  p.fit_ehmos <- out_ehmos$sims.list$p.fit
  p.fitnew_ehmos <- out_ehmos$sims.list$p.fitnew
  bp_ehmos[r] <- sum(p.fit_ehmos > p.fitnew_ehmos)/length(p.fit_ehmos)
  
  # specific
  for (i in 1:N){
    p.fit.sp_ehmos <- out_ehmos$sims.list$p.fit.sp[,i]
    p.fitnew.sp_ehmos <- out_ehmos$sims.list$p.fitnew.sp[,i]
    bp.sp_ehmos[i, r] <- sum(p.fit.sp_ehmos > p.fitnew.sp_ehmos)/length(p.fit.sp_ehmos)
  }
}
apply(bp.sp_ehmos, 1, range)

# for glmm
bp_glmm <- vector("numeric", length = R)
bp.sp_glmm <- matrix(0, nrow = N, ncol = R)
for (r in 1:R){
  # general
  out_glmm <- res.glmm_asym[[r]]
  p.fit_glmm <- out_glmm$sims.list$p.fit
  p.fitnew_glmm <- out_glmm$sims.list$p.fitnew
  bp_glmm[r] <- sum(p.fit_glmm > p.fitnew_glmm)/length(p.fit_glmm)
  
  # specific
  for (i in 1:N){
    p.fit.sp_glmm <- out_glmm$sims.list$p.fit.sp[,i]
    p.fitnew.sp_glmm <- out_glmm$sims.list$p.fitnew.sp[,i]
    bp.sp_glmm[i, r] <- sum(p.fit.sp_glmm > p.fitnew.sp_glmm)/length(p.fit.sp_glmm)
  }
}

apply(bp.sp_glmm, 1, function(x)(sum(x < 0.05 | x > 0.95)))
apply(bp.sp_ehmos, 1, function(x)(sum(x < 0.05 | x > 0.95)))

#### 3. 'Posterior' of SRCs ----
# Take as examples species 1 and 6 (best and worst cases)
sims <- res.ehmos_asym[[1]]$sims.list
x <- X_A1
x.sc <- scale(X_A1)[, 1]
psi.hat1 <- matrix(0, nrow = length(x.sc), ncol = length(sims$opt.mean))
for (m in 1:length(sims$opt.mean)){
  psi.hat1[, m] <- plogis(sims$b0[m, 1, 1] - ((x.sc - sims$opt[m, 1])^2) / (2 * sims$tol[m, 1, 1]^2))
}

psi.hat6 <- matrix(0, nrow = length(x.sc), ncol = length(sims$opt.mean))
for (m in 1:length(sims$opt.mean)){
  psi.hat6[, m] <- plogis(sims$b0[m, 6, 1] - ((x.sc - sims$opt[m, 6])^2) / (2 * sims$tol[m, 6, 1]^2))
}

plot(x.sc, psi.hat1[, 1], type = "l", ylim = c(0, 1))
for (m in 2:length(sims$opt.mean)){
 lines(x.sc, psi.hat1[, m]) 
  lines(x.sc, psi.hat6[, m], col ="red")
}


#### 4. Empiric estimation of occupancy probability to detect asymmetry ----
# test for species 6
cat_alt <- seq(2100, 2600, length.out = 10)
diff_cat <- cat_alt[2] - cat_alt[1]
prop.diff <- matrix(0, nrow = R, ncol = N)
for (i in 1:N){
  for (r in 1:R){
    obs <- data.frame("x" = X_A2,
                      "y" = Z[, i, 1, 2, r]) # first period, so should be skewed on the left of the optimum
    obs %>%
      dplyr::mutate("cat_x" = cut(x, breaks = cat_alt, labels = 1:(length(cat_alt)-1)),
                    "cat_x" = as.numeric(cat_x)) %>%
      dplyr::group_by(cat_x) %>%
      dplyr::summarise_at("y", .funs = function(x)sum(x)/length(x)) -> prop.obs
    if (nrow(prop.obs[prop.obs$y == max(prop.obs$y), "cat_x"]) == 2){
      cat_max <- prop.obs[prop.obs$y == max(prop.obs$y), "cat_x"]
      prop.obs %>%
        dplyr::ungroup() %>%
        dplyr::filter(cat_x != as.numeric(cat_max[1, ])) %>%
        dplyr::mutate("left" = cat_x <= as.numeric(cat_max[1, ])) %>%
        dplyr::group_by(left) %>%
        dplyr::summarise_at("y", sum) -> prop.sum
    }else{
      cat_max <- as.numeric(prop.obs[prop.obs$y == max(prop.obs$y), ]$cat_x) %>% median
      prop.obs %>%
        dplyr::ungroup() %>%
        dplyr::filter(cat_x != cat_max) %>%
        dplyr::mutate("left" = cat_x < cat_max) %>%
        dplyr::group_by(left) %>%
        dplyr::summarise_at("y", sum) -> prop.sum
    }
    
    prop.diff[r, i] <- as.numeric(prop.sum[1, 2] - prop.sum[2, 2])
  }
}



#### density curve of altitudes' occupied
obs %>% dplyr::filter(y == 1) -> pres
plot(density(pres$x))
mtext(print(moments::skewness(pres$x)))

skew <-  matrix(0, nrow = R, ncol = N)
for (i in 1:N){
  for (r in 1:R){
    obs <- data.frame("x" = X_A2,
                      "y" = Z[, i, 1, 2, r]) # first period, so should be skewed on the left of the optimum
    obs %>% dplyr::filter(y == 1) -> pres
    skew[r, i] <- moments::skewness(pres$x)
  }
}
apply(skew > 0.5 | skew < -0.5, 2, sum)

ks.p <-  matrix(0, nrow = R, ncol = N)
for (i in 1:N){
  for (r in 1:R){
    obs <- data.frame("x" = X_A2,
                      "y" = Z[, i, 1, 2, r]) # first period, so should be skewed on the left of the optimum
    obs2 <- data.frame("x" = X_A2,
                      "y" = Z[, i, 2, 2, r])
    obs %>% dplyr::filter(y == 1) -> pres
    obs2 %>% dplyr::filter(y == 1) -> pres2
    ks.p[r, i] <- ks.test((pres$x), (pres2$x))$p.value
  }
}

run.hof <- FALSE
if (run.hof == TRUE){
  shape.hat <-  matrix(0, nrow = R, ncol = N)
  for (i in 1:N){
    for (r in 1:R){
      obs <- data.frame("x" = X_A2,
                        "y" = Z[, i, 1, 2, r]) # first period, so should be skewed on the left of the optimum
      obs2 <- data.frame("x" = X_A2,
                         "y" = Z[, i, 2, 2, r])
      eHOF::HOF(occ = obs$y, grad = obs$x) -> hof
      shape.hat[r, i] <- eHOF::pick.model(hof)
    }
  }
  save(shape.hat, file = here::here("results", "simu_asym_shape_hof.RData"))
}else{
  load(here::here("results", "simu_asym_shape_hof.RData"))
}
  

#### Try to detect asymmetry with GAMs and border estimates (Heegaard 2002) ####
library(mgcv)
diff_bound.gam <-  matrix(0, nrow = R, ncol = N)
for (i in 1:N){
  for (r in 1:R){
    obs <- data.frame("x" = X_A2,
                      "y" = Z[, i, 1, 2, r]) # first period, so should be skewed on the left of the optimum
    gam1 <- gam(y ~ s(x, k = 10), family = binomial(link = "logit"), data = obs)
    d <- data.frame("x" = 1000:3000,
                    "y.hat" = predict.gam(gam1, data.frame("x" = 1000:3000), type = "response"))
    
    # compute out border estimates
    fr <- d[d$y.hat > max(d$y.hat) * exp(-2),]$x
    ob1 <- fr[1]
    ob2 <- fr[length(fr)]
    # diff between out borders and optimum
    opt <- d[d$y.hat == max(d$y.hat), ]$x
    diff_bound.gam[r, i] <- abs(opt - ob1) - abs(opt - ob2)
    # in right skewed distribution we assume negative diff (ob1 closer to opt than ob2)
  }
}
apply(diff_bound.gam, 2, mean) # seems to work in average but not replicate by replicate

## Try to have confidence intervals around the diff between boundaries and optimum
# with bootstrap to then have a test of significance
shape_gam <- function(data, i){ # function used for the bootstrap
  mod <- gam(y ~ s(x, k = 5), family = binomial(link = "logit"), data = data[i, ])
  d <- data.frame("x" = 1000:3000,
                  "y.hat" = predict.gam(mod, data.frame("x" = 1000:3000), type = "response"))
  
  fr <- d[d$y.hat > max(d$y.hat) * exp(-2),]$x
  ob1 <- fr[1]
  ob2 <- fr[length(fr)]
  opt <- d[d$y.hat == max(d$y.hat), ]$x
  return(abs(opt - ob1) - abs(opt - ob2))
}

run_gam <- FALSE
if (run_gam == TRUE){
  gam_shape <-  matrix("0", nrow = R, ncol = N)
  for (i in 1:N){
    for (r in 1:R){
      bootstrap_gam <- boot::boot(data = obs, statistic = shape_gam, R = 100)
      ci <- boot::boot.ci(boot.out=bootstrap_gam,type=c("norm"))
      
      if (ci$normal[2] < 0 & ci$normal[3] > 0){ # if CI include zero we can not
        # conclude about asymmetry
        gam_shape[r, i] <- "IV"
      }else{
        gam_shape[r, i] <- "V"
      }
    }
  }
  gam_shape # always indicating symmetry due to large intervals
  # but could be sensitive to gam and bootstrap specifications
}


