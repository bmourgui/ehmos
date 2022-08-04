                      ######################################%
                      #### Sensitivity analysis results ####%
                      ######################################%

#Reminder: Here I try to assess the sensitivity of the model developed to the choice of prior distributions
source(here::here("analysis", "appli01_formatting-data.R"))
alti.sc <- (alti-mean(alti))/sd(alti)

m.x <- mean(alti)
sd.x <- sd(alti)

load(here::here("results/appli_out-ehmos.RData")) #first run of the application, with usual priors and convergence problem for tau.tol
out_def <- out_appli_ehmos
data <- out$data

load(here::here("results/sensi_out_EHMOS_flat_prior.RData"))
out_flat <- out
load(here::here("results/sensi_out_EHMOS_weak_prior.RData"))
out_weak <- out
load(here::here("results/sensi_out_EHMOS_inf_prior.RData"))
out_inf <- out

#### Saving posterior samples for the 'community-parameters' ####
##Initial model
mu.a1_def <- out_def$sims.list$alpha.mean[,1] #mu.alpha period 1 for the first defel
mu.a2_def <- out_def$sims.list$alpha.mean[,2] #mu.alpha period 2 for the first defel
tau.a_def <- out_def$sims.list$tau.alpha #tau.alpha for the first defel, same for both period
sigma.a_def <- 1/sqrt(tau.a_def) #standard deviation deduced from the precision parameter tau

mu.T1_def <- out_def$sims.list$tol.mean[,1] #mu.T period 1 for the first defel
mu.T2_def <- out_def$sims.list$tol.mean[,2] #mu.T period 2 for the first defel
tau.T_def <- out_def$sims.list$tau.tol #tau.T for the first defel, same for both period
sigma.T_def <- 1/sqrt(tau.T_def) #standard deviation deduced from the precision parameter tau

mu.opt_def <- out_def$sims.list$opt.mean #mu.opt period 1 for the first defel
tau.opt_def <- out_def$sims.list$tau.opt #tau.opt for the first defel, same for both period
sigma.opt_def <- 1/sqrt(tau.opt_def) #standard deviation deduced from the precision parameter tau

mu.shift_def <- out_def$sims.list$shift.mean #mu.shift period 1 for the first defel
tau.shift_def <- out_def$sims.list$tau.shift #tau.shift for the first defel, same for both period
sigma.shift_def <- 1/sqrt(tau.shift_def) #standard deviation deduced from the precision parameter tau

##Initial model
mu.a1_weak <- out_weak$sims.list$alpha.mean[,1] #mu.alpha period 1 for the first weakel
mu.a2_weak <- out_weak$sims.list$alpha.mean[,2] #mu.alpha period 2 for the first weakel
tau.a_weak <- out_weak$sims.list$tau.alpha #tau.alpha for the first weakel, same for both period
sigma.a_weak <- 1/sqrt(tau.a_weak) #standard deviation deduced from the precision parameter tau

mu.T1_weak <- out_weak$sims.list$tol.mean[,1] #mu.T period 1 for the first weakel
mu.T2_weak <- out_weak$sims.list$tol.mean[,2] #mu.T period 2 for the first weakel
tau.T_weak <- out_weak$sims.list$tau.tol #tau.T for the first weakel, same for both period
sigma.T_weak <- 1/sqrt(tau.T_weak) #standard deviation deduced from the precision parameter tau

mu.opt_weak <- out_weak$sims.list$opt.mean #mu.opt period 1 for the first weakel
tau.opt_weak <- out_weak$sims.list$tau.opt #tau.opt for the first weakel, same for both period
sigma.opt_weak <- 1/sqrt(tau.opt_weak) #standard deviation deduced from the precision parameter tau

mu.shift_weak <- out_weak$sims.list$shift.mean #mu.shift period 1 for the first weakel
tau.shift_weak <- out_weak$sims.list$tau.shift #tau.shift for the first weakel, same for both period
sigma.shift_weak <- 1/sqrt(tau.shift_weak) #standard deviation deduced from the precision parameter tau

##second model
mu.a1_flat <- out_flat$sims.list$alpha.mean[,1] #mu.alpha period 1 for the first flatel
mu.a2_flat <- out_flat$sims.list$alpha.mean[,2] #mu.alpha period 2 for the first flatel
tau.a_flat <- out_flat$sims.list$tau.alpha #tau.alpha for the first flatel, same for both period
sigma.a_flat <- 1/sqrt(tau.a_flat) #standard deviation deduced from the precision parameter tau

mu.T1_flat <- out_flat$sims.list$tol.mean[,1] #mu.T period 1 for the first flatel
mu.T2_flat <- out_flat$sims.list$tol.mean[,2] #mu.T period 2 for the first flatel
tau.T_flat <- out_flat$sims.list$tau.tol #tau.T for the first flatel, same for both period
sigma.T_flat <- 1/sqrt(tau.T_flat) #standard deviation deduced from the precision parameter tau

mu.opt_flat <- out_flat$sims.list$opt.mean #mu.opt period 1 for the first flatel
tau.opt_flat <- out_flat$sims.list$tau.opt #tau.opt for the first flatel, same for both period
sigma.opt_flat <- 1/sqrt(tau.opt_flat) #standard deviation deduced from the precision parameter tau

mu.shift_flat <- out_flat$sims.list$shift.mean #mu.shift period 1 for the first flatel
tau.shift_flat <- out_flat$sims.list$tau.shift #tau.shift for the first flatel, same for both period
sigma.shift_flat <- 1/sqrt(tau.shift_flat) #standard deviation deduced from the precision parameter tau


##third model
mu.a1_inf <- out_inf$sims.list$alpha.mean[,1] #mu.alpha period 1 for the first infel
mu.a2_inf <- out_inf$sims.list$alpha.mean[,2] #mu.alpha period 2 for the first infel
tau.a_inf <- out_inf$sims.list$tau.alpha #tau.alpha for the first infel, same for both period
sigma.a_inf <- 1/sqrt(tau.a_inf) #standard deviation deduced from the precision parameter tau

mu.T1_inf <- out_inf$sims.list$tol.mean[,1] #mu.T period 1 for the first infel
mu.T2_inf <- out_inf$sims.list$tol.mean[,2] #mu.T period 2 for the first infel
tau.T_inf <- out_inf$sims.list$tau.tol #tau.T for the first infel, same for both period
sigma.T_inf <- 1/sqrt(tau.T_inf) #standard deviation deduced from the precision parameter tau

mu.opt_inf <- out_inf$sims.list$opt.mean #mu.opt period 1 for the first infel
tau.opt_inf <- out_inf$sims.list$tau.opt #tau.opt for the first infel, same for both period
sigma.opt_inf <- 1/sqrt(tau.opt_inf) #standard deviation deduced from the precision parameter tau

mu.shift_inf <- out_inf$sims.list$shift.mean #mu.shift period 1 for the first infel
tau.shift_inf <- out_inf$sims.list$tau.shift #tau.shift for the first infel, same for both period
sigma.shift_inf <- 1/sqrt(tau.shift_inf) #standard deviation deduced from the precision parameter tau


#### Computing prior distribution used for each community-parameter ####

# Model 1 and 3 #
n <- length(out_weak$sims.list$tau.alpha) #length of posterior sample
mu.a_prior.def <- rnorm(n, mean = 0, sd = 1/sqrt(0.001))
mu.T_prior.def <- rnorm(n, mean = 0, sd = 1/sqrt(0.001))
mu.opt_prior.def <- rnorm(n, mean = 0, sd = 1/sqrt(0.001))
mu.shift_prior.def <- rnorm(n, mean = 0, sd = 1/sqrt(0.001))

tau.prior <- rgamma(n, 0.1, 0.1)

# Model 2 #
mu.a_prior.flat <- rlogis(n, 0, 1)
mu.T_prior.flat <- runif(n, 0, 1.12)
mu.opt_prior.flat <- runif(n, -5.55, 3.45)
mu.shift_prior.flat <- runif(n, -2.8, 2.8)

mu.a_prior.weak <- rnorm(n, mean = -0.6, sd = 1)
mu.T_prior.weak <- msm::rmenorm(n, mean = 0.53, sd = 0.42, lower = 0)
mu.opt_prior.weak <- rnorm(n, mean = 0, sd = 2.8)
mu.shift_prior.weak <- rnorm(n, mean = 0.28, sd = 1.13)

mu.a_prior.inf <- rnorm(n, mean = -0.85, sd = 0.4)
mu.T_prior.inf <- rnorm(n, mean = 0.53, sd = 0.04)
mu.opt_prior.inf <- rnorm(n, mean = -0.48, sd = 0.28)
mu.shift_prior.inf <- rnorm(n, mean = 0.28, sd = 0.05)


#### Graphs comparing posteriors between them and with priors ####
## Alpha1
posterior.def <- density(mu.shift_def)
posterior.flat <- density(mu.shift_flat)
posterior.weak <- density(mu.shift_weak)
#posterior.inf <- density(mu.shift_inf)

prior.def <- density(mu.shift_prior.def)
prior.flat <- density(mu.shift_prior.flat)
prior.weak <- density(mu.shift_prior.weak)
#prior.inf <- density(mu.shift_prior.inf)

plot(x = posterior.def$x * sd.x , 
     y = posterior.def$y, 
     xlim = c(-200, 600),
     ylab = "Scaled density", # all densities are scaled to have max 1
     xlab = "Community shift (m)",
     type ="l", 
     main = "mu_shift")
lines(x = posterior.flat$x * sd.x, y = posterior.flat$y, col="blue")
lines(x = posterior.weak$x * sd.x, y = posterior.weak$y, col="darkgreen")
#lines(x = posterior.inf$x * sd.x, y = posterior.inf$y / max(posterior.inf$y), col="purple")
lines(x = prior.def$x * sd.x, y = prior.def$y, lty = 2)
lines(x = prior.flat$x * sd.x, y = prior.flat$y, col="blue", lty = 2)
lines(x = prior.weak$x * sd.x, y = prior.weak$y, col="darkgreen", lty = 2)
#lines(x = prior.inf$x * sd.x, y = prior.inf$y / max(prior.inf$y), col="purple", lty = 2)


par(mfrow = c(4, 3))
for (i in 1:24){
  
  posterior.def.sp <- density(out_def$sims.list$shift[,i])
  posterior.flat.sp <- density(out_flat$sims.list$shift[,i])
  posterior.weak.sp <- density(out_weak$sims.list$shift[,i])
  posterior.inf.sp <- density(out_inf$sims.list$shift[,i])
  
  
  plot(x = posterior.def.sp$x * sd.x , 
       y = posterior.def.sp$y / max(posterior.def.sp$y), 
       xlim = c(-400, 800),
       ylab = "Scaled density", # all densities are scaled to have max 1
       xlab = "",
       type ="l", 
       main = paste0("Species ", i))
  lines(x = posterior.flat.sp$x * sd.x, y = posterior.flat.sp$y / max(posterior.flat.sp$y), col="blue")
  lines(x = posterior.weak.sp$x * sd.x, y = posterior.weak.sp$y / max(posterior.weak.sp$y), col="darkgreen")
  lines(x = posterior.inf.sp$x * sd.x, y = posterior.inf.sp$y / max(posterior.inf.sp$y), col="purple")
  
  abline(v = out_def$mean$shift[i] * sd.x, lty=2)
  abline(v = out_flat$mean$shift[i] * sd.x, lty=2, col ="blue")
  abline(v = out_weak$mean$shift[i] * sd.x, lty=2, col ="darkgreen")
  abline(v = out_inf$mean$shift[i] * sd.x, lty=2, col ="purple")
  
}

points(out_flat$sims.list$shift[,i], jitter(rep(-0.05,length(out_flat$sims.list$shift[,i])),1.5))
abline(v = out_flat$mean$shift[i], lty=2)
lines(density(out_weak$sims.list$shift[,i]), col="red")
points(out_weak$sims.list$shift[,i], jitter(rep(-0.1,length(out_flat$sims.list$shift[,i])),1.5), col="red")
abline(v = out_weak$mean$shift[i], lty=2, col="red")
lines(density(out_inf$sims.list$shift[,i]), col="green")
points(out_inf$sims.list$shift[,i], jitter(rep(-0.15,length(out_flat$sims.list$shift[,i])),1.5), col="green")
abline(v = out_inf$mean$shift[i], lty=2, col="green")