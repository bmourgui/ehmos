                      ######################################%
                      #### Sensitivity analysis results ####%
                      ######################################%

#Reminder: Here I try to assess the sensitivity of the model developed to the choice of prior distributions

load("~/Bureau/PhD/Orthoptera/CC/Data/out_application.RData") #first run of the application, with usual priors and convergence problem for tau.tol
load("~/Bureau/PhD/Orthoptera/CC/Analysis/out_appli_jamil.RData")
out <- out_appli_jamil_sensi1
load("~/Bureau/PhD/Orthoptera/CC/Analysis/out_appli_jamil_sensi1.RData") #second run with usual priors, no convergence issue
load("~/Bureau/PhD/Orthoptera/CC/Analysis/out_appli_jamil_sensi2.RData") #run with my own priors, 'more informative' but with the same form (normal for mu's and gamma for tau's)
load("~/Bureau/PhD/Orthoptera/CC/Analysis/out_appli_jamil_sensi3.RData")

#saving outputs within objects with smaller names
out1 <- out_appli_jamil_sensi1
out2 <- out_appli_jamil_sensi2
out3 <- out_appli_jamil_sensi3

#### Saving posterior samples for the 'community-parameters' ####
##Initial model
mu.a1_mod <- out$sims.list$b0.mean[,1] #mu.alpha period 1 for the first model
mu.a2_mod <- out$sims.list$b0.mean[,2] #mu.alpha period 2 for the first model
tau.a_mod <- out$sims.list$tau.b0 #tau.alpha for the first model, same for both period
sigma.a_mod <- 1/sqrt(tau.a_mod) #standard deviation deduced from the precision parameter tau

mu.T1_mod <- out$sims.list$tol.mean[,1] #mu.T period 1 for the first model
mu.T2_mod <- out$sims.list$tol.mean[,2] #mu.T period 2 for the first model
tau.T_mod <- out$sims.list$tau.tol #tau.T for the first model, same for both period
sigma.T_mod <- 1/sqrt(tau.T_mod) #standard deviation deduced from the precision parameter tau

mu.opt_mod <- out$sims.list$opt.mean #mu.opt period 1 for the first model
tau.opt_mod <- out$sims.list$tau.opt #tau.opt for the first model, same for both period
sigma.opt_mod <- 1/sqrt(tau.opt_mod) #standard deviation deduced from the precision parameter tau

mu.shift_mod <- out$sims.list$shift.mean #mu.shift period 1 for the first model
tau.shift_mod <- out$sims.list$tau.shift #tau.shift for the first model, same for both period
sigma.shift_mod <- 1/sqrt(tau.shift_mod) #standard deviation deduced from the precision parameter tau

##Initial model
mu.a1_mod1 <- out1$sims.list$alpha.mean[,1] #mu.alpha period 1 for the first mod1el
mu.a2_mod1 <- out1$sims.list$alpha.mean[,2] #mu.alpha period 2 for the first mod1el
tau.a_mod1 <- out1$sims.list$tau.alpha #tau.alpha for the first mod1el, same for both period
sigma.a_mod1 <- 1/sqrt(tau.a_mod1) #standard deviation deduced from the precision parameter tau

mu.T1_mod1 <- out1$sims.list$tol.mean[,1] #mu.T period 1 for the first mod1el
mu.T2_mod1 <- out1$sims.list$tol.mean[,2] #mu.T period 2 for the first mod1el
tau.T_mod1 <- out1$sims.list$tau.tol #tau.T for the first mod1el, same for both period
sigma.T_mod1 <- 1/sqrt(tau.T_mod1) #standard deviation deduced from the precision parameter tau

mu.opt_mod1 <- out1$sims.list$opt.mean #mu.opt period 1 for the first mod1el
tau.opt_mod1 <- out1$sims.list$tau.opt #tau.opt for the first mod1el, same for both period
sigma.opt_mod1 <- 1/sqrt(tau.opt_mod1) #standard deviation deduced from the precision parameter tau

mu.shift_mod1 <- out1$sims.list$shift.mean #mu.shift period 1 for the first mod1el
tau.shift_mod1 <- out1$sims.list$tau.shift #tau.shift for the first mod1el, same for both period
sigma.shift_mod1 <- 1/sqrt(tau.shift_mod1) #standard deviation deduced from the precision parameter tau

##second model
mu.a1_mod2 <- out2$sims.list$alpha.mean[,1] #mu.alpha period 1 for the second model
mu.a2_mod2 <- out2$sims.list$alpha.mean[,2] #mu.alpha period 2 for the second model
tau.a_mod2 <- out2$sims.list$tau.alpha #tau.alpha for the second model, same for both period
sigma.a_mod2 <- 1/sqrt(tau.a_mod2) #standard deviation deduced from the precision parameter tau

mu.T1_mod2 <- out2$sims.list$tol.mean[,1] #mu.T period 1 for the second model
mu.T2_mod2 <- out2$sims.list$tol.mean[,2] #mu.T period 2 for the second model
tau.T_mod2 <- out2$sims.list$tau.tol #tau.T for the second model, same for both period
sigma.T_mod2 <- 1/sqrt(tau.T_mod2) #standard deviation deduced from the precision parameter tau

mu.opt_mod2 <- out2$sims.list$opt.mean #mu.opt period 1 for the second model
tau.opt_mod2 <- out2$sims.list$tau.opt #tau.opt for the second model, same for both period
sigma.opt_mod2 <- 1/sqrt(tau.opt_mod2) #standard deviation deduced from the precision parameter tau

mu.shift_mod2 <- out2$sims.list$shift.mean #mu.shift period 1 for the second model
tau.shift_mod2 <- out2$sims.list$tau.shift #tau.shift for the second model, same for both period
sigma.shift_mod2 <- 1/sqrt(tau.shift_mod2) #standard deviation deduced from the precision parameter tau


##third model
mu.a1_mod3 <- out3$sims.list$alpha.mean[,1] #mu.alpha period 1 for the third model
mu.a2_mod3 <- out3$sims.list$alpha.mean[,2] #mu.alpha period 2 for the second model
tau.a_mod3 <- out3$sims.list$tau.alpha #tau.alpha for the third model, same for both period
sigma.a_mod3 <- 1/sqrt(tau.a_mod3) #standard deviation deduced from the precision parameter tau

mu.T1_mod3 <- out3$sims.list$tol.mean[,1] #mu.T period 1 for the third model
mu.T2_mod3 <- out3$sims.list$tol.mean[,2] #mu.T period 2 for the second model
tau.T_mod3 <- out3$sims.list$tau.tol #tau.T for the third model, same for both period
sigma.T_mod3 <- 1/sqrt(tau.T_mod3) #standard deviation deduced from the precision parameter tau

mu.opt_mod3 <- out3$sims.list$opt.mean #mu.opt period 1 for the third model
tau.opt_mod3 <- out3$sims.list$tau.opt #tau.opt for the third model, same for both period
sigma.opt_mod3 <- 1/sqrt(tau.opt_mod3) #standard deviation deduced from the precision parameter tau

mu.shift_mod3 <- out3$sims.list$shift.mean #mu.shift period 1 for the third model
tau.shift_mod3 <- out3$sims.list$tau.shift #tau.shift for the third model, same for both period
sigma.shift_mod3 <- 1/sqrt(tau.shift_mod3) #standard deviation deduced from the precision parameter tau


#### Computing prior distribution used for each community-parameter ####

# Model 1 and 3 #
n <- length(out2$sims.list$tau.b0) #length of posterior sample
mu.a_prior <- rnorm(n, mean = 0, sd = 1/sqrt(0.001))
mu.T_prior <- rnorm(n, mean = 0, sd = 1/sqrt(0.001))
mu.opt_prior <- rnorm(n, mean = 0, sd = 1/sqrt(0.001))
mu.shift_prior <- rnorm(n, mean = 0, sd = 1/sqrt(0.001))

tau.a_prior <- rgamma(n, 0.1, 0.1)
sigma.a_prio <- 1/sqrt(tau.a_prior)
tau.T_prior <- rgamma(n, 0.1, 0.1)
sigma.T_prio <- 1/sqrt(tau.T_prior)
tau.opt_prior <- rgamma(n, 0.1, 0.1)
sigma.opt_prio <- 1/sqrt(tau.opt_prior)
tau.shift_prior <- rgamma(n, 0.1, 0.1)
sigma.shift_prio <- 1/sqrt(tau.shift_prior)

# Model 2 #
mu.a_prior2 <- rnorm(n, mean = 0.5, sd = 1/sqrt(0.0625))
mu.T_prior2 <- rnorm(n, mean = 0.25, sd = 1)
mu.opt_prior2 <- rnorm(n, mean = -1, sd = 1/sqrt(0.002))
mu.shift_prior2 <- rnorm(n, mean = 0.2, sd = 1/sqrt(0.01))

tau.a_prior2 <- rgamma(n, 0.1, 0.1)
sigma.a_prio <- 1/sqrt(tau.a_prior2)
tau.T_prior2 <- rgamma(n, 0.1, 0.1)
sigma.T_prio <- 1/sqrt(tau.T_prior2)
tau.opt_prior2 <- rgamma(n, 0.1, 0.1)
sigma.opt_prio <- 1/sqrt(tau.opt_prior2)
tau.shift_prior2 <- rgamma(n, 0.1, 0.1)
sigma.shift_prio <- 1/sqrt(tau.shift_prior2)


#### Graphs comparing posteriors between them and with priors ####
## Alpha1
plot(density(mu.a1_mod1), main = "Alpha1")
lines(density(mu.a1_mod2), col="darkgreen")
lines(density(mu.a1_mod3), col="purple")
lines(density(mu.a1_mod), col="blue")
par(new=T)
plot(density(mu.a_prior), lty=2, 
     ylim=c(min(density(mu.a_prior)$y), max(density(mu.a_prior)$y)), 
     xlim=c(min(density(mu.a1_mod1)$y), max(density(mu.a1_mod1)$y)),
     axes=F, main="", ylab="", xlab="")
par(new=T)
plot(density(mu.a_prior2), lty=2, 
     ylim=c(min(density(mu.a_prior2)$y), max(density(mu.a_prior2)$y)), 
     xlim=c(min(density(mu.a1_mod1)$y), max(density(mu.a1_mod1)$y)),
     axes=F, main="", ylab="", xlab="")

## Alpha2
plot(density(mu.a2_mod), main = "Alpha2")
lines(density(mu.a2_mod2), col="darkgreen")
lines(density(mu.a2_mod3), col="purple")
par(new=T)
plot(density(mu.a_prior), lty=2, 
     ylim=c(min(density(mu.a_prior)$y), max(density(mu.a_prior)$y)), 
     xlim=c(min(density(mu.a2_mod1)$y), max(density(mu.a2_mod1)$y)),
     axes=F, main="", ylab="", xlab="")
par(new=T)
plot(density(mu.a_prior2), lty=2, 
     ylim=c(min(density(mu.a_prior2)$y), max(density(mu.a_prior2)$y)), 
     xlim=c(min(density(mu.a1_mod1)$y), max(density(mu.a1_mod1)$y)),
     axes=F, main="", ylab="", xlab="")

## Tol1
plot(density(mu.T1_mod1), main = "Tol1", 
     ylim=c(min(c(density(mu.T1_mod1)$y, density(mu.T1_mod2)$y, density(mu.T1_mod3)$y)), 
            max(c(density(mu.T1_mod1)$y, density(mu.T1_mod2)$y, density(mu.T1_mod3)$y))))#,
     #xlim=c(min(density(mu.T_prior2)$x), max(density(mu.T1_mod2)$x)))
lines(density(mu.T1_mod2), col="darkgreen")
lines(density(mu.T1_mod3), col="purple")
lines(density(mu.T1_mod), col="red")
par(new=T)
plot(density(mu.T_prior), lty=2, 
     ylim=c(min(density(mu.T_prior)$y), max(density(mu.T_prior)$y)), 
     axes=F, main="", ylab="", xlab="",
     xlim=c(min(density(mu.T_prior2)$x), max(density(mu.T1_mod2)$x)))
par(new=T)
plot(density(mu.T_prior2), lty=2, 
     ylim=c(min(density(mu.T_prior2)$y), max(density(mu.T_prior2)$y)),
     axes=F, main="", ylab="", xlab="", col="purple",
     xlim=c(min(density(mu.T_prior2)$x), max(density(mu.T1_mod2)$x)))

## Tol2
plot(density(mu.T2_mod1), main = "Tol2", 
     ylim=c(min(c(density(mu.T2_mod1)$y, density(mu.T2_mod2)$y, density(mu.T2_mod3)$y)), 
            max(c(density(mu.T2_mod1)$y, density(mu.T2_mod2)$y, density(mu.T2_mod3)$y))))
lines(density(mu.T2_mod2), col="darkgreen")
lines(density(mu.T2_mod3), col="purple")
par(new=T)
plot(density(mu.T_prior), lty=2, 
     ylim=c(min(density(mu.T_prior)$y), max(density(mu.T_prior)$y)), 
     xlim=c(min(density(mu.T2_mod1)$y), max(density(mu.T2_mod1)$y)),
     axes=F, main="", ylab="", xlab="")


## Opt
plot(density(mu.opt_mod1), main = "opt", 
     ylim=c(min(c(density(mu.opt_mod1)$y, density(mu.opt_mod2)$y, density(mu.opt_mod3)$y)), 
            max(c(density(mu.opt_mod1)$y, density(mu.opt_mod2)$y, density(mu.opt_mod3)$y))))
lines(density(mu.opt_mod2), col="darkgreen")
lines(density(mu.opt_mod3), col="purple")
par(new=T)
plot(density(mu.T_prior), lty=2, 
     ylim=c(min(density(mu.T_prior)$y), max(density(mu.T_prior)$y)), 
     xlim=c(min(density(mu.opt_mod1)$y), max(density(mu.opt_mod1)$y)),
     axes=F, main="", ylab="", xlab="")


## Shift
plot(density(mu.shift_mod1), main = "shift", 
     ylim=c(min(c(density(mu.shift_mod1)$y, density(mu.shift_mod2)$y, density(mu.shift_mod3)$y)), 
            max(c(density(mu.shift_mod1)$y, density(mu.shift_mod2)$y, density(mu.shift_mod3)$y))))
lines(density(mu.shift_mod2), col="darkgreen")
lines(density(mu.shift_mod3), col="purple")


plot(density(out1$sims.list$b0.mean[,1]), ylim=c(-0.15,max(max(density(out1$sims.list$b0.mean[,1])$y),max(density(out2$sims.list$b0.mean[,1])$y))))
points(out1$sims.list$b0.mean[,1], jitter(rep(-0.05,length(out1$sims.list$b0.mean[,1])),1.5))
lines(density(out2$sims.list$b0.mean[,1]), col="red")
points(out2$sims.list$b0.mean[,1], jitter(rep(-0.1,length(out1$sims.list$b0.mean[,1])),1.5), col="red")
lines(density(out$sims.list$b0.mean[,1]), col="green")
points(out$sims.list$b0.mean[,1], jitter(rep(-0.15,length(out1$sims.list$b0.mean[,1])),1.5), col="green")


plot(density(out1$sims.list$b0.mean[,2]), ylim=c(-0.15,max(max(density(out1$sims.list$b0.mean[,2])$y),max(density(out2$sims.list$b0.mean[,2])$y))))
points(out1$sims.list$b0.mean[,2], jitter(rep(-0.05,length(out1$sims.list$b0.mean[,2])),1.5))
lines(density(out2$sims.list$b0.mean[,2]), col="red")
points(out2$sims.list$b0.mean[,2], jitter(rep(-0.1,length(out1$sims.list$b0.mean[,2])),1.5), col="red")
lines(density(out$sims.list$b0.mean[,2]), col="green")
points(out$sims.list$b0.mean[,2], jitter(rep(-0.15,length(out1$sims.list$b0.mean[,2])),1.5), col="green")


plot(density(out1$sims.list$opt.mean), ylim=c(-0.15,max(max(density(out1$sims.list$opt.mean)$y),max(density(out2$sims.list$opt.mean)$y))))
points(out1$sims.list$opt.mean, jitter(rep(-0.05,length(out1$sims.list$opt.mean)),1.5))
lines(density(out2$sims.list$opt.mean), col="red")
points(out2$sims.list$opt.mean, jitter(rep(-0.1,length(out1$sims.list$opt.mean)),1.5), col="red")
lines(density(out$sims.list$opt.mean), col="green")
points(out$sims.list$opt.mean, jitter(rep(-0.15,length(out1$sims.list$opt.mean)),1.5), col="green")


plot(density(out1$sims.list$shift.mean), ylim=c(-0.15,max(max(density(out1$sims.list$shift.mean)$y),max(density(out2$sims.list$shift.mean)$y))))
points(out1$sims.list$shift.mean, jitter(rep(-0.05,length(out1$sims.list$shift.mean)),1.5))
lines(density(out2$sims.list$shift.mean), col="red")
points(out2$sims.list$shift.mean, jitter(rep(-0.1,length(out1$sims.list$shift.mean)),1.5), col="red")
lines(density(out$sims.list$shift.mean), col="green")
points(out$sims.list$shift.mean, jitter(rep(-0.15,length(out1$sims.list$shift.mean)),1.5), col="green")


plot(density(out1$sims.list$tol.mean[,1]), ylim=c(-0.15,max(max(density(out1$sims.list$tol.mean[,1])$y),max(density(out2$sims.list$tol.mean[,1])$y))))
points(out1$sims.list$tol.mean[,1], jitter(rep(-0.05,length(out1$sims.list$tol.mean[,1])),1.5))
lines(density(out2$sims.list$tol.mean[,1]), col="red")
points(out2$sims.list$tol.mean[,1], jitter(rep(-0.1,length(out1$sims.list$tol.mean[,1])),1.5), col="red")
lines(density(out$sims.list$tol.mean[,1]), col="green")
points(out$sims.list$tol.mean[,1], jitter(rep(-0.15,length(out1$sims.list$tol.mean[,1])),1.5), col="green")


plot(density(out1$sims.list$tol.mean[,2]), ylim=c(-0.15,max(max(density(out1$sims.list$tol.mean[,2])$y),max(density(out2$sims.list$tol.mean[,2])$y))))
points(out1$sims.list$tol.mean[,2], jitter(rep(-0.05,length(out1$sims.list$tol.mean[,2])),1.5))
lines(density(out2$sims.list$tol.mean[,2]), col="red")
points(out2$sims.list$tol.mean[,2], jitter(rep(-0.1,length(out1$sims.list$tol.mean[,2])),1.5), col="red")
lines(density(out$sims.list$tol.mean[,2]), col="green")
points(out$sims.list$tol.mean[,2], jitter(rep(-0.15,length(out1$sims.list$tol.mean[,2])),1.5), col="green")



plot(density(out1$sims.list$tau.tol), ylim=c(-0.015,max(max(density(c(out$samples[[1]][,"tau.tol"],out$samples[[2]][,"tau.tol"]))$y),max(density(out2$sims.list$tau.tol)$y))))
points(out1$sims.list$tau.tol, jitter(rep(-0.005,length(out1$sims.list$tau.tol)),1.5))
lines(density(out2$sims.list$tau.tol), col="red")
points(out2$sims.list$tau.tol, jitter(rep(-0.01,length(out1$sims.list$tau.tol)),1.5), col="red")
lines(density(c(out$samples[[2]][,"tau.tol"],out$samples[[3]][,"tau.tol"])), col="green")
points(c(out$samples[[1]][,"tau.tol"],out$samples[[2]][,"tau.tol"]), jitter(rep(-0.015,length(c(out$samples[[1]][,"tau.tol"],out$samples[[2]][,"tau.tol"]))),1.5), col="green")

plot(density(out1$sims.list$tau.opt), ylim=c(-0.15,max(max(density(out1$sims.list$tau.opt)$y),max(density(out2$sims.list$tau.opt)$y))))
points(out1$sims.list$tau.opt, jitter(rep(-0.05,length(out1$sims.list$tau.opt)),1.5))
lines(density(out2$sims.list$tau.opt), col="red")
points(out2$sims.list$tau.opt, jitter(rep(-0.1,length(out1$sims.list$tau.opt)),1.5), col="red")
lines(density(out$sims.list$tau.opt), col="green")
points(out$sims.list$tau.opt, jitter(rep(-0.15,length(out1$sims.list$tau.opt)),1.5), col="green")


plot(density(out1$sims.list$tau.shift), ylim=c(-0.15,max(max(density(out1$sims.list$tau.shift)$y),max(density(out2$sims.list$tau.shift)$y))))
points(out1$sims.list$tau.shift, jitter(rep(-0.05,length(out1$sims.list$tau.shift)),1.5))
lines(density(out2$sims.list$tau.shift), col="red")
points(out2$sims.list$tau.shift, jitter(rep(-0.1,length(out1$sims.list$tau.shift)),1.5), col="red")
lines(density(out$sims.list$tau.shift), col="green")
points(out$sims.list$tau.shift, jitter(rep(-0.15,length(out1$sims.list$tau.shift)),1.5), col="green")

plot(density(out1$sims.list$tau.b0), ylim=c(-0.15,max(max(density(out1$sims.list$tau.b0)$y),max(density(out2$sims.list$tau.b0)$y))))
points(out1$sims.list$tau.b0, jitter(rep(-0.05,length(out1$sims.list$tau.b0)),1.5))
lines(density(out2$sims.list$tau.b0), col="red")
points(out2$sims.list$tau.b0, jitter(rep(-0.1,length(out1$sims.list$tau.b0)),1.5), col="red")
lines(density(out$sims.list$tau.b0), col="green")
points(out$sims.list$tau.b0, jitter(rep(-0.15,length(out1$sims.list$tau.b0)),1.5), col="green")



for (i in 1:24){
  plot(density(out1$sims.list$b0[,i]), ylim=c(-0.15,max(max(density(out1$sims.list$b0[,i])$y),max(density(out2$sims.list$b0[,i])$y))))
  points(out1$sims.list$b0[,i], jitter(rep(-0.05,length(out1$sims.list$b0[,i])),1.5))
  abline(v = out1$mean$b0[i], lty=2)
  lines(density(out2$sims.list$b0[,i]), col="red")
  points(out2$sims.list$b0[,i], jitter(rep(-0.1,length(out1$sims.list$b0[,i])),1.5), col="red")
  abline(v = out2$mean$b0[i], lty=2, col="red")
  lines(density(out$sims.list$b0[,i]), col="green")
  points(out$sims.list$b0[,i], jitter(rep(-0.15,length(out1$sims.list$b0[,i])),1.5), col="green")
  abline(v = out$mean$b0[i], lty=2, col="green")
}

for (i in 1:24){
  plot(density(out1$sims.list$opt[,i]), ylim=c(-0.15,max(max(density(out1$sims.list$opt[,i])$y),max(density(out2$sims.list$opt[,i])$y))))
  points(out1$sims.list$opt[,i], jitter(rep(-0.05,length(out1$sims.list$opt[,i])),1.5))
  abline(v = out1$mean$opt[i], lty=2)
  lines(density(out2$sims.list$opt[,i]), col="red")
  points(out2$sims.list$opt[,i], jitter(rep(-0.1,length(out1$sims.list$opt[,i])),1.5), col="red")
  abline(v = out2$mean$opt[i], lty=2, col="red")
  lines(density(out$sims.list$opt[,i]), col="green")
  points(out$sims.list$opt[,i], jitter(rep(-0.15,length(out1$sims.list$opt[,i])),1.5), col="green")
  abline(v = out$mean$opt[i], lty=2, col="green")
}

for (i in 1:24){
  plot(density(out1$sims.list$tol[,i,1]), ylim=c(-0.15,max(max(density(out1$sims.list$tol[,i,1])$y),max(density(out2$sims.list$tol[,i,1])$y))))
  points(out1$sims.list$tol[,i,1], jitter(rep(-0.05,length(out1$sims.list$tol[,i,1])),1.5))
  abline(v = out1$mean$tol[i], lty=2)
  lines(density(out2$sims.list$tol[,i,1]), col="red")
  points(out2$sims.list$tol[,i,1], jitter(rep(-0.1,length(out1$sims.list$tol[,i,1])),1.5), col="red")
  abline(v = out2$mean$tol[i], lty=2, col="red")
  lines(density(out$sims.list$tol[,i,1]), col="green")
  points(out$sims.list$tol[,i,1], jitter(rep(-0.15,length(out1$sims.list$tol[,i,1])),1.5), col="green")
  abline(v = out$mean$tol[i], lty=2, col="green")
}

for (i in 1:24){
  plot(density(out1$sims.list$b0[,i]), ylim=c(-0.15,max(max(density(out1$sims.list$b0[,i])$y),max(density(out2$sims.list$b0[,i])$y))))
  points(out1$sims.list$b0[,i], jitter(rep(-0.05,length(out1$sims.list$b0[,i])),1.5))
  abline(v = out1$mean$b0[i], lty=2)
  lines(density(out2$sims.list$b0[,i]), col="red")
  points(out2$sims.list$b0[,i], jitter(rep(-0.1,length(out1$sims.list$b0[,i])),1.5), col="red")
  abline(v = out2$mean$b0[i], lty=2, col="red")
  lines(density(out$sims.list$b0[,i]), col="green")
  points(out$sims.list$b0[,i], jitter(rep(-0.15,length(out1$sims.list$b0[,i])),1.5), col="green")
  abline(v = out$mean$b0[i], lty=2, col="green")
}