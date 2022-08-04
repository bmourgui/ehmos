library(magrittr)
library(jagsUI)

#########################%
#### Data importation ####
#########################%

a <- c(-5, 0, 5)
o <- c(0, 2000, 3200) # highest peak at 3,143 m
o.sc <- (o - m.x) / sd.x
w <- c(0, 900, 3000)
t <- 0.5 * w / sqrt(2*(a - log(0.05/0.95)))
t[1] <- 0
t.sc <- t / sd.x
d <- c(-1000, 200, 1000)
d.sc <- d / sd.x
x <- -1000:5000
x.sc <- (x - m.x) / sd.x

psi.inf <- plogis(a[1] - ((x - o[1])^2)/(2 * t[1] ^ 2))
psi.mean <- plogis(a[2] - ((x - o[2])^2)/(2 * (t[2] ^ 2)))
psi.sup <- plogis(a[3] - ((x - o[3])^2)/(2 * (t[3] ^ 2)))

psi.sc.inf <- plogis(a[1] - ((x.sc - o.sc[1])^2)/(2 * t.sc[1] ^ 2))
psi.sc.mean <- plogis(a[2] - ((x.sc - o.sc[2])^2)/(2 * (t.sc[2] ^ 2)))
psi.sc.sup <- plogis(a[3] - ((x.sc - o.sc[3])^2)/(2 * (t.sc[3] ^ 2)))

plot(x, psi.mean, type = "l", ylim = c(0, 1))
lines(x, psi.inf)
lines(x, psi.sup)
lines(x, psi.sc.mean, col ="red")
lines(x, psi.sc.inf, col ="red")
lines(x, psi.sc.sup, col ="red")

source(here::here("analysis", "appli01_formatting-data.R"))
alti.sc <- (alti-mean(alti))/sd(alti)

m.x <- mean(alti)
sd.x <- sd(alti)


##### Visualization of prior distributions ####
par(mfrow = c(2, 2))
## Optimum
default.prior <- rnorm(100000, 0, 1/sqrt(0.001))
default.prior.unscaled <- default.prior * sd.x + m.x
dens.default.prior <- density(default.prior.unscaled)

flat.prior <- runif(100000, min = (0 - m.x) / sd.x, max = (3200 - m.x) / sd.x)
flat.prior.unscaled <- flat.prior * sd.x + m.x
dens.flat.prior <- density(flat.prior.unscaled)

weak.prior <- rnorm(100000, mean = (2000 - m.x) / sd.x, sd = 1000 / sd.x)
weak.prior.unscaled <- weak.prior * sd.x + m.x
dens.weak.prior <- density(weak.prior.unscaled)

plot(x = dens.weak.prior$x, y = dens.weak.prior$y, 
     xlim = c(-2000, 6000), 
     xlab = "Elevation",
     ylab = "density",
     type = "l", 
     main = expression(theta))
lines(x = dens.flat.prior$x, y = dens.flat.prior$y, col = "blue")
lines(x = dens.default.prior$x, y = dens.default.prior$y, col = "red")


## Shift
default.prior <- rnorm(100000, 0, 1/sqrt(0.001))
default.prior.unscaled <- default.prior * sd.x
dens.default.prior <- density(default.prior.unscaled)

flat.prior <- runif(100000, min = (- 1000 / sd.x), max = 1000 / sd.x)
flat.prior.unscaled <- flat.prior * sd.x
dens.flat.prior <- density(flat.prior.unscaled)

weak.prior <- rnorm(100000, mean = 100 / sd.x, sd = 400 / sd.x)
weak.prior.unscaled <- weak.prior * sd.x
dens.weak.prior <- density(weak.prior.unscaled)

plot(x = dens.weak.prior$x, y = dens.weak.prior$y, 
     xlim = c(-2000, 2000), 
     xlab = "Elevation",
     ylab = "",
     type = "l", 
     main = expression(delta))
lines(x = dens.flat.prior$x, y = dens.flat.prior$y, col = "blue")
lines(x = dens.default.prior$x, y = dens.default.prior$y, col = "red")

## Tolerance
default.prior <- msm::rmenorm(100000, 0, 1/sqrt(0.001), lower = 0)
default.prior.unscaled <- default.prior * sd.x
dens.default.prior <- density(default.prior.unscaled)

flat.prior <- runif(100000, min = 0, max = 400 / sd.x)
flat.prior.unscaled <- flat.prior * sd.x
dens.flat.prior <- density(flat.prior.unscaled)

weak.prior <- msm::rmenorm(100000, mean = 190 / sd.x, sd = 150 / sd.x, lower = 0)
weak.prior.unscaled <- weak.prior * sd.x
dens.weak.prior <- density(weak.prior.unscaled)

plot(x = dens.weak.prior$x, y = dens.weak.prior$y, 
     xlim = c(0, 1000), 
     xlab = "Elevation",
     ylab = "density",
     type = "l", 
     main = expression(tau))
lines(x = dens.flat.prior$x, y = dens.flat.prior$y, col = "blue")
lines(x = dens.default.prior$x, y = dens.default.prior$y, col = "red")


## Maximum probability
default.prior <- rnorm(100000, 0, 1/sqrt(0.001))
default.prior.unscaled <- plogis(default.prior)
dens.default.prior <- density(default.prior.unscaled)

flat.prior <- rlogis(100000, 0, 1) # dlogis in JAGS, advised by @northrup2018
flat.prior.unscaled <- plogis(flat.prior)
dens.flat.prior <- density(flat.prior.unscaled)

weak.prior <- rnorm(100000, mean = -0.6, sd = 1) # 0.35 on probability scale
weak.prior.unscaled <- plogis(weak.prior)
dens.weak.prior <- density(weak.prior.unscaled)

plot(x = dens.default.prior$x, y = dens.default.prior$y, 
     xlim = c(0, 1), 
     xlab = "Probability",
     ylab = "",
     type = "l", 
     col = "red",
     main = expression(psi[max]))
lines(x = dens.flat.prior$x, y = dens.flat.prior$y, col = "blue")
lines(x = dens.weak.prior$x, y = dens.weak.prior$y)


###########################################################%
#### Hyper parameters selection based on assumptions ####
###########################################################%

## Optimum parameter
# we assumed a community mean around 1800 m
# and a community sd around 400 m in ecological scale.

# We simulate a normal distribution with those parameters
# to find the parameters at model scale
mu.opt <- rnorm(10000, 1800, 400)
mu.opt.sc <- (mu.opt - m.x) / sd.x
mean(mu.opt.sc) ; 1/(sd(mu.opt.sc)^2)

# so we will set prior distributions centered in -0.48 for the mean
# and in 0.8 for the precision parameter

# priors for the community mean, normal distributions centered
# on 1800 m (at ecological scale), with increasing precision
# and then transformed on model scale
prior.opt1 <- rnorm(10000, 1800, 1000)
prior.opt2 <- rnorm(10000, 1800, 500)
prior.opt3 <- rnorm(10000, 1800, 100)

sc.prior.opt1 <- (prior.opt1 - m.x) / sd.x
sc.prior.opt2 <- (prior.opt2 - m.x) / sd.x
sc.prior.opt3 <- (prior.opt3 - m.x) / sd.x

# obtain precision
mean(sc.prior.opt1) ; mean(sc.prior.opt2) ; mean(sc.prior.opt3)
1/(sd(sc.prior.opt1) ^ 2) ; 1/(sd(sc.prior.opt2) ^ 2) ; 1/(sd(sc.prior.opt3) ^ 2)

# priors for the community variance, gamma distributions, centered
# on 0.8 (at model scale) with increasing precision


## Tolerance parameter
# First make assumption on width
# 
width <- rnorm(100000, 1000, 300)
tol <- 0.5 * width / sqrt(2 * (0.5 - log(0.05/0.95)))
mean(tol) ; sd(tol) ; 1/(sd(tol) ^ 2)

prior.tol1 <- rnorm(10000, 190, 100)
prior.tol2 <- rnorm(10000, 190, 50)
prior.tol3 <- rnorm(10000, 190, 10)

sc.prior.tol1 <- (prior.tol1) / sd.x
sc.prior.tol2 <- (prior.tol2) / sd.x
sc.prior.tol3 <- (prior.tol3 ) / sd.x

# obtain precision
1/(sd(sc.prior.tol1) ^ 2) ; 1/(sd(sc.prior.tol2) ^ 2) ; 1/(sd(sc.prior.tol3) ^ 2)


## Shift parameter
# We could assumed mostly upward shifts, and shifts ranged between -400 +800 m
shift <- rnorm(10000, 200, 200)
sc.shift <- shift / sd.x
mean(sc.shift) ; sd(sc.shift)


#### Visualize "prior predictive check" ####
## I simulate SRC based on 100 draw of prior distributions
it <- 100

m.opt.flat <- runif(it, -5.55, 3.45)
m.opt.weak <- rnorm(it, mean = 0, sd = 2.8)
m.opt.default <- rnorm(it, mean = 0, sd = 1/sqrt(0.001))

m.shift.flat <- runif(it, -2.28, 2.28)
m.shift.weak <- rnorm(it, mean = 0.28, sd = 1.13)
m.shift.default <- rnorm(it, mean = 0, sd = 1/sqrt(0.001))

m.tol.flat <- runif(it, -0, 1.12)
m.tol.weak <- msm::rmenorm(it, mean = 0.53, sd = 0.42, lower = 0)
m.tol.default <- rnorm(it, mean = 0, sd = 1/sqrt(0.001))

m.alpha.flat <- rlogis(it, 0, 1)
m.alpha.weak <- rnorm(it, mean = -0.6, sd = 1)
m.alpha.default <- rnorm(it, mean = 0, sd = 1/sqrt(0.001))

x <- seq(-10000, 10000, 1)
sc.x <- (x - m.x)/sd.x

m.psi.flat <- matrix(0, nrow = length(x), ncol = length(m.alpha.flat))
m.psi.weak <- matrix(0, nrow = length(x), ncol = length(m.alpha.flat))
m.psi.default <- matrix(0, nrow = length(x), ncol = length(m.alpha.flat))

for (j in 1:length(x)){
  m.psi.flat[j, ] <- plogis(m.alpha.flat - 0.5 * ((sc.x[j] - m.opt.flat)^2) / (m.tol.flat ^ 2))
  m.psi.weak[j, ] <- plogis(m.alpha.weak - 0.5 * ((sc.x[j] - m.opt.weak)^2) / (m.tol.weak ^ 2))
  m.psi.default[j, ] <- plogis(m.alpha.default - 0.5 * ((sc.x[j] - m.opt.default)^2) / (m.tol.default ^ 2))
}

par(mfrow = c(3, 1))
plot(x = x, 
     y = m.psi.flat[, 1], 
     type = "l", 
     col = "darkgreen", 
     ylim = c(0, 1), 
     xlim = c(-1000, 5000),
     main = "Flat priors",
     ylab = expression(psi))
for (i in 1:100){
  lines(x = x, y = m.psi.flat[, i], col = "darkgreen")
}
plot(x = x, 
     y = m.psi.weak[, 1], 
     type = "l", 
     col = "blue", 
     ylim = c(0, 1), 
     xlim = c(-1000, 5000),
     main = "Weak priors",
     ylab = expression(psi))
for (i in 1:100){
  lines(x = x, y = m.psi.weak[, i], col = "blue")
}
plot(x = x, 
     y = m.psi.default[, 1], 
     type = "l", 
     col = "red", 
     ylim = c(0, 1), 
     xlim = c(-1000, 5000),
     main = "Default priors",
     ylab = expression(psi))
for (i in 1:100){
  lines(x = x, y = m.psi.default[, i], col = "red")
}


#### Model with informative priors ####

cat("
    model{

  for (i in 1:n){
    for (t in 1:S){
      for (j in 1:J){
        
        z[j,i,t] ~ dbern(psi.bound[j,i,t])
        psi.bound[j,i,t] <- max(0,min(1,psi[j,i,t]))
        logit(psi[j,i,t]) <- alpha[i,t] - 0.5*pow((alti[j]-(opt[i]+shift[i]*survey[t]))/(tol[i,t]),2) + a.site[j]
        
      }
    }
  }

for (i in 1:n){

    opt[i] ~ dnorm(opt.mean, tau.opt)
    shift[i] ~ dnorm(shift.mean, tau.shift)

    for (t in 1:S){

    	tol[i,t] ~ dnorm(tol.mean[t], tau.tol)T(0.01,)
    	alpha[i,t] ~ dnorm(alpha.mean[t], tau.alpha)

    }
  }


for (t in 1:S){
    tol.mean[t] ~ dunif(0, 1.12)
    alpha.mean[t] ~ dlogis(0, 1)
  }

    opt.mean ~ dunif(-5.55, 3.45)
    shift.mean ~ dunif(-2.8, 2.8)


  for (j in 1:J){
    a.site[j] ~ dnorm(0, tau.a)
  }

  tau.opt ~ dgamma(0.1,0.1)
  tau.tol ~ dgamma(0.1,0.1)
  tau.alpha ~ dgamma(0.1,0.1)
  tau.shift ~ dgamma(0.1,0.1)
  tau.a ~ dgamma(0.1,0.1)

}", file="mod_EHMOS_flat_prior.txt")

#### Model with flatormative priors ####

cat("
    model{

  for (i in 1:n){
    for (t in 1:S){
      for (j in 1:J){
        
        z[j,i,t] ~ dbern(psi.bound[j,i,t])
        psi.bound[j,i,t] <- max(0,min(1,psi[j,i,t]))
        logit(psi[j,i,t]) <- alpha[i,t] - 0.5*pow((alti[j]-(opt[i]+shift[i]*survey[t]))/(tol[i,t]),2) + a.site[j]
        
      }
    }
  }

for (i in 1:n){

    opt[i] ~ dnorm(opt.mean, tau.opt)
    shift[i] ~ dnorm(shift.mean, tau.shift)

    for (t in 1:S){

    	tol[i,t] ~ dnorm(tol.mean[t], tau.tol)T(0.01,)
    	alpha[i,t] ~ dnorm(alpha.mean[t], tau.alpha)

    }
  }


for (t in 1:S){
    tol.mean[t] ~ dnorm(0.53, 1/(0.42 ^ 2))T(0.01,)
    alpha.mean[t] ~ dnorm(-0.6, 1)
  }

    opt.mean ~ dnorm(0, 1/(2.8^2))
    shift.mean ~ dnorm(0.28, 1/(1.13^2))


  for (j in 1:J){
    a.site[j] ~ dnorm(0, tau.a)
  }

  tau.opt ~ dgamma(0.1,0.1)
  tau.tol ~ dgamma(0.1,0.1)
  tau.alpha ~ dgamma(0.1,0.1)
  tau.shift ~ dgamma(0.1,0.1)
  tau.a ~ dgamma(0.1,0.1)

}", file="mod_EHMOS_weak_prior.txt")


#### Specify the data ####
sp.data = list(n=n, J=J, S=S, 
               z=Z, 
               alti=alti.sc, survey=c(0,1))

#### Specify the parameters to be monitored ####
params_ehmos = c('alpha.mean','tau.alpha','alpha',
               'opt.mean','tau.opt', 'opt', 
               'tol.mean','tau.tol', 'tol', 
               'shift.mean','tau.shift', 'shift',
               'tau.a','a.site'
)

inits_ehmos = function() {
  psi.meanGuess=runif(1, 0, 1)
  list(a.site=rnorm(J), 
       pmax.mean=runif(2, 0, 1),
       opt.mean=rnorm(1, 0, 0.1),
       shift.mean=rnorm(1, 0, 0.1),
       tol.mean=runif(2,0.1,0.6))
}


library(foreach)
cl <- parallel::makeForkCluster(2)
doParallel::registerDoParallel(cl)

foreach (s=1:2) %dopar% {
    set.seed(64)
    prior <- c("EHMOS_flat_prior", "EHMOS_weak_prior")
    out <- autojags(data = sp.data,
                    inits = inits_ehmos,
                    parameters.to.save = params_ehmos,
                    model.file = paste0("mod_", prior[s], ".txt"),
                    n.chains = 3,
                    n.adapt = 5000,
                    n.thin = 15,
                    n.burnin = 15000,
                    iter.increment = 15000,
                    DIC = TRUE,
                    store.data = TRUE, #keep initial values
                    parallel = F,
                    max.iter = 250000)
    save(out, file = here::here("results", paste0("sensi_out_", prior[s], ".RData"))) # in case of crash
  }


load(here::here("results/appli_out-ehmos.RData")) #first run of the application, with usual priors and convergence problem for tau.tol
out_def <- out_appli_ehmos
data <- out_def$data

load(here::here("results/sensi_out_EHMOS_flat_prior.RData"))
out_flat <- out
load(here::here("results/sensi_out_EHMOS_weak_prior.RData"))
out_weak <- out

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

posterior.opt_def <- density(mu.opt_def)
posterior.opt_flat <- density(mu.opt_flat)
posterior.opt_weak <- density(mu.opt_weak)

posterior.shift_def <- density(mu.shift_def)
posterior.shift_flat <- density(mu.shift_flat)
posterior.shift_weak <- density(mu.shift_weak)

posterior.tol1_def <- density(mu.T1_def)
posterior.tol1_flat <- density(mu.T1_flat)
posterior.tol1_weak <- density(mu.T1_weak)

posterior.tol2_def <- density(mu.T2_def)
posterior.tol2_flat <- density(mu.T2_flat)
posterior.tol2_weak <- density(mu.T2_weak)

posterior.a1_def <- density(mu.a1_def)
posterior.a1_flat <- density(mu.a1_flat)
posterior.a1_weak <- density(mu.a1_weak)

posterior.a2_def <- density(mu.a2_def)
posterior.a2_flat <- density(mu.a2_flat)
posterior.a2_weak <- density(mu.a2_weak)

par(mfrow = c(2, 3))
plot(x = posterior.opt_def$x * sd.x + m.x, 
     y = posterior.opt_def$y, 
     xlim = c(min(posterior.opt_def$x * sd.x + m.x), max(posterior.opt_def$x * sd.x + m.x)),
     ylim = c(0, max(posterior.opt_weak$y)),
     ylab = "Density", # all densities are scaled to have max 1
     xlab = "Elevation (m)",
     type ="l", 
     main = expression(mu[theta]))
abline(v = median(posterior.opt_def$x * sd.x + m.x), lty = 2)
lines(x = posterior.opt_flat$x * sd.x + m.x, y = posterior.opt_flat$y, col="blue")
abline(v = median(posterior.opt_flat$x * sd.x + m.x), lty = 2, col="blue")
lines(x = posterior.opt_weak$x * sd.x + m.x, y = posterior.opt_weak$y, col="darkgreen")
abline(v = median(posterior.opt_weak$x * sd.x + m.x), lty = 2, col="darkgreen")

plot(x = posterior.tol1_def$x * sd.x, 
     y = posterior.tol1_def$y, 
     xlim = c(min(posterior.tol1_def$x * sd.x), max(posterior.tol1_def$x * sd.x)),
     ylim = c(0, max(posterior.tol1_weak$y)),
     ylab = "Density", # all densities are scaled to have max 1
     xlab = "Elevation (m)",
     type ="l", 
     main = expression(mu[tau[1]]))
abline(v = median(posterior.tol1_def$x * sd.x), lty = 2)
lines(x = posterior.tol1_flat$x * sd.x, y = posterior.tol1_flat$y, col="blue")
abline(v = median(posterior.tol1_flat$x * sd.x), lty = 2, col="blue")
lines(x = posterior.tol1_weak$x * sd.x, y = posterior.tol1_weak$y, col="darkgreen")
abline(v = median(posterior.tol1_weak$x * sd.x), lty = 2, col="darkgreen")

plot(x = plogis(posterior.a1_def$x), 
     y = posterior.a1_def$y, 
     xlim = c(0, 1),
     ylim = c(0, max(posterior.a1_weak$y)),
     ylab = "Density", # all densities are scaled to have max 1
     xlab = "Maximum probability",
     type ="l", 
     main = expression(mu[alpha[1]]))
abline(v = median(plogis(posterior.a1_def$x)), lty = 2)
lines(x = plogis(posterior.a1_flat$x), y = posterior.a1_flat$y, col="blue")
abline(v = median(plogis(posterior.a1_flat$x)), lty = 2, col="blue")
lines(x = plogis(posterior.a1_weak$x), y = posterior.a1_weak$y, col="darkgreen")
abline(v = median(plogis(posterior.a1_weak$x)), lty = 2, col="darkgreen")

plot(x = posterior.shift_def$x * sd.x, 
     y = posterior.shift_def$y, 
     xlim = c(min(posterior.shift_def$x * sd.x), max(posterior.shift_def$x * sd.x)),
     ylim = c(0, max(posterior.shift_weak$y)),
     ylab = "Density", # all densities are scaled to have max 1
     xlab = "Elevation (m)",
     type ="l", 
     main = expression(mu[delta]))
abline(v = median(posterior.shift_def$x * sd.x), lty = 2)
lines(x = posterior.shift_flat$x * sd.x, y = posterior.shift_flat$y, col="blue")
abline(v = median(posterior.shift_flat$x * sd.x), lty = 2, col="blue")
lines(x = posterior.shift_weak$x * sd.x, y = posterior.shift_weak$y, col="darkgreen")
abline(v = median(posterior.shift_weak$x * sd.x), lty = 2, col="darkgreen")

plot(x = posterior.tol2_def$x * sd.x, 
     y = posterior.tol2_def$y, 
     xlim = c(min(posterior.tol2_def$x * sd.x), max(posterior.tol2_def$x * sd.x)),
     ylim = c(0, max(posterior.tol2_weak$y)),
     ylab = "Density", # all densities are scaled to have max 1
     xlab = "Elevation (m)",
     type ="l", 
     main = expression(mu[tau[2]]))
abline(v = median(posterior.tol2_def$x * sd.x), lty = 2)
lines(x = posterior.tol2_flat$x * sd.x, y = posterior.tol2_flat$y, col="blue")
abline(v = median(posterior.tol2_flat$x * sd.x), lty = 2, col="blue")
lines(x = posterior.tol2_weak$x * sd.x, y = posterior.tol2_weak$y, col="darkgreen")
abline(v = median(posterior.tol2_weak$x * sd.x), lty = 2, col="darkgreen")

plot(x = plogis(posterior.a2_def$x), 
     y = posterior.a2_def$y, 
     xlim = c(0, 1),
     ylim = c(0, max(posterior.a2_weak$y)),
     ylab = "Density", # all densities are scaled to have max 1
     xlab = "Maximum probability",
     type ="l", 
     main = expression(mu[alpha[2]]))
abline(v = median(plogis(posterior.a2_def$x)), lty = 2)
lines(x = plogis(posterior.a2_flat$x), y = posterior.a2_flat$y, col="blue")
abline(v = median(plogis(posterior.a2_flat$x)), lty = 2, col="blue")
lines(x = plogis(posterior.a2_weak$x), y = posterior.a2_weak$y, col="darkgreen")
abline(v = median(plogis(posterior.a2_weak$x)), lty = 2, col="darkgreen")



par(mfrow = c(4, 3))
for (i in 1:24){
  
  posterior.def.sp <- density(out_def$sims.list$shift[,i])
  posterior.flat.sp <- density(out_flat$sims.list$shift[,i])
  posterior.weak.sp <- density(out_weak$sims.list$shift[,i])
  
  full.sp.name <- as.character(dimnames(data$Z)[[2]][i])
  split.sp.name <- strsplit(full.sp.name, split = " ")[[1]]
  sp.name <- paste(substr(split.sp.name[1], 1, 2), split.sp.name[2], sep = ". ")
  
  plot(x = posterior.def.sp$x * sd.x , 
       y = posterior.def.sp$y, 
       xlim = c(-400, 800),
       ylim = c(0, max(posterior.weak.sp$y)),
       ylab = "Density", # all densities are scaled to have max 1
       xlab = "",
       type ="l", 
       main = sp.name)
  lines(x = posterior.flat.sp$x * sd.x, y = posterior.flat.sp$y , col="blue")
  lines(x = posterior.weak.sp$x * sd.x, y = posterior.weak.sp$y , col="darkgreen")
  
  abline(v = out_def$mean$shift[i] * sd.x, lty=2)
  abline(v = out_flat$mean$shift[i] * sd.x, lty=2, col ="blue")
  abline(v = out_weak$mean$shift[i] * sd.x, lty=2, col ="darkgreen")
  
}


par(mfrow = c(4, 3))
for (i in 1:24){
  
  posterior.def.sp <- density(out_def$sims.list$opt[,i])
  posterior.flat.sp <- density(out_flat$sims.list$opt[,i])
  posterior.weak.sp <- density(out_weak$sims.list$opt[,i])
  
  full.sp.name <- as.character(dimnames(data$Z)[[2]][i])
  split.sp.name <- strsplit(full.sp.name, split = " ")[[1]]
  sp.name <- paste(substr(split.sp.name[1], 1, 2), split.sp.name[2], sep = ". ")
  
  plot(x = posterior.def.sp$x * sd.x + m.x , 
       y = posterior.def.sp$y, 
       xlim = c(0, 3200),
       ylim = c(0, max(posterior.weak.sp$y)),
       ylab = "Density", # all densities are scaled to have max 1
       xlab = "",
       type ="l", 
       main = sp.name)
  lines(x = posterior.flat.sp$x * sd.x + m.x, y = posterior.flat.sp$y , col="blue")
  lines(x = posterior.weak.sp$x * sd.x + m.x, y = posterior.weak.sp$y , col="darkgreen")
  
  abline(v = out_def$mean$opt[i] * sd.x + m.x, lty=2)
  abline(v = out_flat$mean$opt[i] * sd.x + m.x, lty=2, col ="blue")
  abline(v = out_weak$mean$opt[i] * sd.x + m.x, lty=2, col ="darkgreen")
  
}


par(mfrow = c(4, 3))
for (i in 1:24){
  
  posterior.def.sp <- density(out_def$sims.list$tol[,i,1])
  posterior.flat.sp <- density(out_flat$sims.list$tol[,i,1])
  posterior.weak.sp <- density(out_weak$sims.list$tol[,i,1])
  
  full.sp.name <- as.character(dimnames(data$Z)[[2]][i])
  split.sp.name <- strsplit(full.sp.name, split = " ")[[1]]
  sp.name <- paste(substr(split.sp.name[1], 1, 2), split.sp.name[2], sep = ". ")
  
  plot(x = posterior.def.sp$x * sd.x , 
       y = posterior.def.sp$y, 
       xlim = c(100, 600),
       ylim = c(0, max(posterior.weak.sp$y)),
       ylab = "Density", # all densities are scaled to have max 1
       xlab = "",
       type ="l", 
       main = sp.name)
  lines(x = posterior.flat.sp$x * sd.x, y = posterior.flat.sp$y , col="blue")
  lines(x = posterior.weak.sp$x * sd.x, y = posterior.weak.sp$y , col="darkgreen")
  
  abline(v = out_def$mean$tol[i,1] * sd.x, lty=2)
  abline(v = out_flat$mean$tol[i,1] * sd.x, lty=2, col ="blue")
  abline(v = out_weak$mean$tol[i,1] * sd.x, lty=2, col ="darkgreen")
  
}


par(mfrow = c(4, 3))
for (i in 1:24){
  
  posterior.def.sp <- density(out_def$sims.list$tol[,i,2])
  posterior.flat.sp <- density(out_flat$sims.list$tol[,i,2])
  posterior.weak.sp <- density(out_weak$sims.list$tol[,i,2])
  
  full.sp.name <- as.character(dimnames(data$Z)[[2]][i])
  split.sp.name <- strsplit(full.sp.name, split = " ")[[1]]
  sp.name <- paste(substr(split.sp.name[1], 1, 2), split.sp.name[2], sep = ". ")
  
  plot(x = posterior.def.sp$x * sd.x , 
       y = posterior.def.sp$y, 
       xlim = c(100, 600),
       ylim = c(0, max(posterior.weak.sp$y)),
       ylab = "Density", # all densities are scaled to have max 1
       xlab = "",
       type ="l", 
       main = sp.name)
  lines(x = posterior.flat.sp$x * sd.x, y = posterior.flat.sp$y , col="blue")
  lines(x = posterior.weak.sp$x * sd.x, y = posterior.weak.sp$y , col="darkgreen")
  
  abline(v = out_def$mean$tol[i,2] * sd.x, lty=2)
  abline(v = out_flat$mean$tol[i,2] * sd.x, lty=2, col ="blue")
  abline(v = out_weak$mean$tol[i,2] * sd.x, lty=2, col ="darkgreen")
  
}


par(mfrow = c(4, 3))
for (i in 1:24){
  
  posterior.def.sp <- density(out_def$sims.list$alpha[,i,1])
  posterior.flat.sp <- density(out_flat$sims.list$alpha[,i,1])
  posterior.weak.sp <- density(out_weak$sims.list$alpha[,i,1])
  
  full.sp.name <- as.character(dimnames(data$Z)[[2]][i])
  split.sp.name <- strsplit(full.sp.name, split = " ")[[1]]
  sp.name <- paste(substr(split.sp.name[1], 1, 2), split.sp.name[2], sep = ". ")
  
  plot(x = plogis(posterior.def.sp$x)  , 
       y = posterior.def.sp$y, 
       xlim = c(0, 1),
       ylim = c(0, max(posterior.weak.sp$y)),
       ylab = "Density", # all densities are scaled to have max 1
       xlab = "",
       type ="l", 
       main = sp.name)
  lines(x = plogis(posterior.flat.sp$x) , y = posterior.flat.sp$y , col="blue")
  lines(x = plogis(posterior.weak.sp$x) , y = posterior.weak.sp$y , col="darkgreen")
  
  abline(v = plogis(out_def$mean$alpha[i,1]) , lty=2)
  abline(v = plogis(out_flat$mean$alpha[i,1]) , lty=2, col ="blue")
  abline(v = plogis(out_weak$mean$alpha[i,1]) , lty=2, col ="darkgreen")
  
}

par(mfrow = c(4, 3))
for (i in 1:24){
  
  posterior.def.sp <- density(out_def$sims.list$alpha[,i,2])
  posterior.flat.sp <- density(out_flat$sims.list$alpha[,i,2])
  posterior.weak.sp <- density(out_weak$sims.list$alpha[,i,2])
  
  full.sp.name <- as.character(dimnames(data$Z)[[2]][i])
  split.sp.name <- strsplit(full.sp.name, split = " ")[[1]]
  sp.name <- paste(substr(split.sp.name[1], 1, 2), split.sp.name[2], sep = ". ")
  
  plot(x = plogis(posterior.def.sp$x)  , 
       y = posterior.def.sp$y, 
       xlim = c(0, 1),
       ylim = c(0, max(posterior.weak.sp$y)),
       ylab = "Density", # all densities are scaled to have max 1
       xlab = "",
       type ="l", 
       main = sp.name)
  lines(x = plogis(posterior.flat.sp$x) , y = posterior.flat.sp$y , col="blue")
  lines(x = plogis(posterior.weak.sp$x) , y = posterior.weak.sp$y , col="darkgreen")
  
  abline(v = plogis(out_def$mean$alpha[i,2]) , lty=2)
  abline(v = plogis(out_flat$mean$alpha[i,2]) , lty=2, col ="blue")
  abline(v = plogis(out_weak$mean$alpha[i,2]) , lty=2, col ="darkgreen")
  
}
