library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(reshape)
library(reshape2)
library(Hmisc)
library(jagsUI)

mean.method2 <- function(z, alti){
  z[z==0] <- NA
  
  sp.alti <- z*alti
  
  test <- t.test(sp.alti[,2], sp.alti[,1])
  
  res <- c("shift"=test$estimate[[1]]-test$estimate[[2]], "lwr"=test$conf.int[1], "median"=NA,"upr"=test$conf.int[2])
}

#########################%
#### Data importation ####
#########################%

source("~/Bureau/PhD/Orthoptera/CC/Analysis/Application_analysis/Code/Data_formatting.R")
alti.sc <- (alti-mean(alti))/sd(alti)

###########################################;
#### Model following Jamil et al., 2014 ####
###########################################;
{

cat("
    model{
  
  for (i in 1:n){
    for (t in 1:S){
      for (j in 1:J){
        Z[j,i,t] ~ dbern(psi.bound[j,i,t])
        psi.bound[j,i,t] <- max(0,min(1,psi[j,i,t]))
        logit(psi[j,i,t]) <- alpha[i,t] - 0.5*pow((alti[j]-(opt[i]+shift[i]*survey[t]))/(tol[i,t]),2) + a.site[j]
      
      #creat simulated dataset to calculate the Bayesain p-value
        Znew[j,i,t] ~ dbern(psi.bound[j,i,t])
        d[j,i,t] <- abs(Z[j,i,t] - psi.bound[j,i,t])
        dnew[j,i,t] <- abs(Znew[j,i,t] - psi.bound[j,i,t])
        d2[j,i,t] <- pow(d[j,i,t],2)
        dnew2[j,i,t] <- pow(dnew[j,i,t],2)
      }
    }
  }
  
  p.fit <- sum(d2[1:J,1:n,1:S])
  p.fitnew <- sum(dnew2[1:J,1:n,1:S])

for (i in 1:n){

    opt[i] ~ dnorm(opt.mean, tau.opt)#T(-4.42, 3.2)
    shift[i] ~ dnorm(shift.mean, tau.shift)#T(-1, 1)

    for (t in 1:S){

    	tol[i,t] ~ dnorm(tol.mean[t], tau.tol)T(0.05,)  
    	alpha[i,t] ~ dnorm(alpha.mean[t], tau.alpha)

    }
  }


for (t in 1:S){
    tol.mean[t] ~ dnorm(0, 0.001)                        
    alpha.mean[t] ~ dnorm(0, 0.001) 
  }
  
    opt.mean ~ dnorm(0, 0.001)
    shift.mean ~ dnorm(0, 0.001)

  
  for (j in 1:J){
    a.site[j] ~ dnorm(0, tau.a)
  }
  
  tau.opt ~ dgamma(0.1,0.1)
  tau.tol ~ dgamma(0.1,0.1)
  tau.alpha ~ dgamma(0.1,0.1)
  tau.shift ~ dgamma(0.1,0.1)
  tau.a ~ dgamma(0.1,0.1)
  
}", file="mod_jamil_shift-glmm_multi.txt")


#### Specify the data ####
sp.data = list(n=n, J=J, S=S, 
               Z=Z, 
               alti=alti.sc, survey=c(0,1))

#### Specify the parameters to be monitored ####
sp.params = c('alpha.mean','tau.alpha','alpha',
              'opt.mean','tau.opt', 'opt', 
              'tol.mean','tau.tol', 'tol', 
              'shift.mean','tau.shift', 'shift',
              'tau.a','a.site',
              'p.fit','p.fitnew'
)

inits_jamil = function() {
  psi.meanGuess=runif(1, 0, 1)
  list(a.site=rnorm(J), 
       alpha.mean=rnorm(2, mean=2, sd=1),
       opt.mean=rnorm(1,1800,300),
       shift.mean=rnorm(1, 150,100),
       tol.mean=rnorm(2,330,50))
}


out_appli_jamil <- autojags(data = sp.data,
                            inits = inits_jamil,
                            parameters.to.save = sp.params,
                            model.file = "mod_jamil_shift-glmm_multi.txt",
                            n.chains = 3,
                            n.adapt = 1000,
                            n.thin = 10,
                            n.burnin = 15000,
                            iter.increment = 15000,
                            DIC = T,
                            store.data = T,
                            Rhat.limit = 1.1,
                            parallel = T,
                            max.iter = 150000)

save(out_appli_jamil, file = "out_appli_jamil.RData")

}




cat("
    model{
    
    #Define prior distributions for community-level model parameters
    b0.mean ~ dunif(0,1)                          #species level effect on occupancy probability
    mu.b0 <- log(b0.mean) - log(1-b0.mean)        #species level effect on occupancy probability on the logit scale
    
    mub1 ~ dnorm(0, 0.001)                        #detection covariable 1 effect on detectability
    mub2 ~ dnorm(0, 0.001)
    mub3 ~ dnorm(0, 0.001)
    mub4 ~ dnorm(0, 0.001)
    mub5 ~ dnorm(0, 0.001)
    
    tau.a ~ dgamma(0.5,0.005)                      #variability of species level effect on detectability    
    tau.b0 ~ dgamma(0.1, 0.1)
    tau.b1 ~ dgamma(0.1,0.1)                      #variability of detection covariable 1 effect on detectability 
    tau.b2 ~ dgamma(0.1,0.1)
    tau.b3 ~ dgamma(0.1,0.1)
    tau.b4 ~ dgamma(0.1,0.1)
    tau.b5 ~ dgamma(0.1,0.1)
    
    for(j in 1:J){
      a[j] ~ dnorm(0, tau.a) #random site effects
    }
    
    for (i in 1:n) {
    
      #Create priors for species i from the community level prior distributions
      b0[i] ~ dnorm(mu.b0, tau.b0) 
      b1[i] ~ dnorm(mub1, tau.b1)               #detection covariable 1 effect on detectability of the species i
      b2[i] ~ dnorm(mub2, tau.b2)
      b3[i] ~ dnorm(mub3, tau.b3)
      b4[i] ~ dnorm(mub4, tau.b4)
      b5[i] ~ dnorm(mub5, tau.b5)
      
      #Create a loop to estimate the Z matrix for species i 
      #at point j.      
        for (t in 1:S) {
          for (j in 1:J){
            Z[j,i,t] ~ dbern(psi.bound[j,i,t])
            psi.bound[j,i,t]<-max(0,min(1,psi[j,i,t]))
            logit(psi[j,i,t]) <- b0[i] + b1[i]*alti[j] + b2[i]*pow(alti[j],2) +
                                 b3[i]*survey[t] + b4[i]*alti[j]*survey[t] + b5[i]*pow(alti[j],2)*survey[t] +
        		                     a[j]
            #creat simulated dataset to calculate the Bayesain p-value
        Znew[j,i,t] ~ dbern(psi.bound[j,i,t])
        d[j,i,t] <- abs(Z[j,i,t] - psi.bound[j,i,t])
        dnew[j,i,t] <- abs(Znew[j,i,t] - psi.bound[j,i,t])
        d2[j,i,t] <- pow(d[j,i,t],2)
        dnew2[j,i,t] <- pow(dnew[j,i,t],2)
      }
    }
  }
  
  p.fit <- sum(d2[1:J,1:n,1:S])
  p.fitnew <- sum(dnew2[1:J,1:n,1:S])
    }
    
", file="mod_glmm_multi.txt")


params_glmm = c('mu.b0','mub1','mub2','mub3','mub4','mub5',
                'tau.b0','tau.b1','tau.b2','tau.b3','tau.b4','tau.b5',
                'tau.a', 
                'b0','b1','b2','b3','b4','b5',
                'a',
                'p.fit','p.fitnew'
)

#### Specify initial values and run the model 
inits_glmm = function() {
  psi.meanGuess=runif(1, 0, 1)
  list(a=rnorm(J), 
       b0.mean=runif(1,0,1),
       mub1=rnorm(1), 
       mub2=rnorm(1), 
       mub3=rnorm(1),
       mub4=rnorm(1),
       mub5=rnorm(1))
}


out_appli_glmm <- autojags(data = sp.data, 
                  inits = inits_glmm, 
                  parameters.to.save = params_glmm, 
                  model.file = "mod_glmm_multi.txt", 
                  n.chains = 3,
                  n.adapt = 1000,
                  n.thin = 10,
                  n.burnin = 15000,
                  iter.increment = 15000,
                  DIC = T,
                  store.data = T,
                  Rhat.limit = 1.1,
                  parallel = T,
                  max.iter = 150000)

save(out_appli_glmm, out_appli_jamil, file = "~/Bureau/PhD/Orthoptera/CC/Analysis/Application_analysis/Results_data/out_models_appli.RData")



#########################################################################%
#### Sensitivity analysis of 'shift-GLMM' with alti in the true scale ####
#########################################################################%

{
  ###################################%
  #### Hyper parameters selection ####
  ###################################%
  
  ## Check for species data range ##
  dataF[,-c(1,2)]*altitudes$altitude.x -> occu_alti
  occu_alti[occu_alti == 0] <- NA
  sp_range <- apply(occu_alti,2,function(x)max(x,na.rm=T))-apply(occu_alti,2,function(x)min(x,na.rm=T))
  sd(sp_range) ; summary(sp_range) ; hist(sp_range)
  # -> # in the data, species range (width) looks like a normal distribution, centered in 1100m, with a sd of 240m
  t1 <- rnorm(10000, 330/sd(altitudes$altitude.x), 70/sd(altitudes$altitude.x))
  z1 <- t1*sqrt(2*(1-log(0.01/0.99)))
  t2 <- rnorm(10000, 330, 350)
  z2 <- t2*sqrt(2*(1-log(0.01/0.99)))
  t3 <- rnorm(10000, 330, 700)
  z3 <- t3*sqrt(2*(1-log(0.01/0.99)))
  plot(density(z1*sd(altitudes$altitude.x)), xlim=c(0,3000), col="blue")
  lines(density(z2), col="green")
  lines(density(z3), col="red")
  lines(density(sp_range))
  # -> # hence tolerance prior distribution for set 1 (informative set) is tau ~ N(330/sd.site, (70/sd.site)^2)
  # -> # set 2 is tau ~ N(330/sd.site, ((70/sd.site)*5)^2)
  # -> # set 3 is tau ~ N(330/sd.site, ((70/sd.site)*10)^2)
  
  ## Check for species optimum ##
  sp_optD <- apply(occu_alti[dataF$campaign=="historic",], 2, function(x)mean(x, na.rm=T))
  sd(sp_optD) ; summary(sp_optD) ; hist(sp_optD)
  # -> # in the data, species optimum average during period 1 is around 1800, and sd around 240
  opt1 <- rnorm(1000, 1800, 240)
  opt2 <- rnorm(1000, 1800, 240*5)
  opt3 <- rnorm(1000, 1800, 240*10)
  plot(density(sp_optD), xlim=c(0,4000))
  lines(density(opt1), col="blue")
  lines(density(opt2), col="green")
  lines(density(opt3), col="red")
  # -> # hence optimum prior distribution for set 1 (informative set) is tau ~ N((1800-m.site)/sd.site, (240-m.site)/sd.site^2)
  # -> # set 2 is tau ~ N((1800-m.site)/sd.site, ((240-m.site)/sd.site*5)^2)
  # -> # set 3 is tau ~ N((1800-m.site)/sd.site, ((240-m.site)/sd.site*10)^2)
  
  ## Check for species shift ##
  sp_shiftD <- as.numeric(apply(Z,2,function(x)mean.method2(x, altitudes$altitude.x))[1,])
  sd(sp_shiftD) ; summary(sp_shiftD) ; hist(sp_shiftD)
  # -> # in the data, species shiftimum average during period 1 is around 120, and sd around 100
  shift1 <- rnorm(1000, 120, 100)
  shift2 <- rnorm(1000, 120, 100*5)
  shift3 <- rnorm(1000, 120, 100*10)
  plot(density(sp_shiftD), xlim=c(-800,800))
  lines(density(shift1*sd(altitudes$altitude.x)), col="blue")
  lines(density(shift2), col="green")
  lines(density(shift3), col="red")
  # -> # hence tolerance prior distribution for set 1 (informative set) is tau ~ N(120/sd.site, (100/sd.site)^2)
  # -> # set 2 is tau ~ N(120/sd.site, ((100/sd.site)*5)^2)
  # -> # set 3 is tau ~ N(120/sd.site, ((100/sd.site)*10)^2)
  
  ## Check for alpha ##
  # no information directly in the data
  # I used informations from the other article, where orthoptera distribution is modelled in the PNM
  # -> # average of maximum probability is around 0.45, with sd = 0.27
  pmax1 <- rnorm(1000, 0.43, 0.27)
  pmax2 <- rnorm(1000, 0.43, 0.27*5)
  pmax3 <- rnorm(1000, 0.43, 0.27*10)
  alpha1 <- rnorm(1000, log(0.43/0.57), -log(0.27/0.73))
  alpha2 <- rnorm(1000, log(0.43/0.57), 2.5)#-log(0.27/0.73)*5)
  alpha3 <- rnorm(1000, log(0.43/0.57), -log(0.27/0.73)*10)
  plot(density(plogis(alpha1)), col="blue", xlim=c(0,1))
  lines(density(plogis(alpha2)), col="green")
  lines(density(plogis(alpha3)), col="red")
  lines(density(pmax1), col="blue", lty=2)
  lines(density(pmax2), col="green", lty=2)
  lines(density(pmax3), col="red", lty=2)
  
  
  
  
  #### Model sensitivity 2 : set1 ####
  
  cat("
    model{
  
  for (i in 1:n){
    for (t in 1:S){
      for (j in 1:J){
        Z[j,i,t] ~ dbern(psi.bound[j,i,t])
        psi.bound[j,i,t] <- max(0,min(1,psi[j,i,t]))
        logit(psi[j,i,t]) <- alpha[i,t] - 0.5*pow((alti[j]-(opt[i]+shift[i]*survey[t]))/(tol[i,t]),2) + a.site[j]
      
      #creat simulated dataset to calculate the Bayesain p-value
        Znew[j,i,t] ~ dbern(psi.bound[j,i,t])
        d[j,i,t] <- abs(Z[j,i,t] - psi.bound[j,i,t])
        dnew[j,i,t] <- abs(Znew[j,i,t] - psi.bound[j,i,t])
        d2[j,i,t] <- pow(d[j,i,t],2)
        dnew2[j,i,t] <- pow(dnew[j,i,t],2)
      }
    }
  }
  
  p.fit <- sum(d2[1:J,1:n,1:S])
  p.fitnew <- sum(dnew2[1:J,1:n,1:S])

for (i in 1:n){

    opt[i] ~ dnorm(opt.mean, tau.opt)#T(-4.42, 3.2)
    shift[i] ~ dnorm(shift.mean, tau.shift)#T(-1, 1)

    for (t in 1:S){

    	tol[i,t] ~ dnorm(tol.mean[t], tau.tol)T(0.05,)  
    	alpha[i,t] ~ dnorm(alpha.mean[t], tau.alpha)

    }
  }


for (t in 1:S){
    tol.mean[t] ~ dnorm(330, 1/70^2) #seems to be way more informative than sensi 1 but in fact the precision parameter (inverse of variance) is still wide for a parameter with the scale of the tolerance                        
    pmax.mean[t] ~ dnorm(0.45,1/0.27^2)T(0,1)
    alpha.mean[t] <- log(pmax.mean[t])-log(1-pmax.mean[t])
  }
  
    opt.mean ~ dnorm(1800, 1/240^2) #we can assume that the mean is more likely to be in lower elevations
    shift.mean ~ dnorm(120, 1/100^2) #we can assume that mean shit is upward, and we could be pretty sure about it

  
  for (j in 1:J){
    a.site[j] ~ dnorm(0, tau.a)
  }
  
  tau.opt ~ dgamma(0.1,0.1)
  tau.tol ~ dgamma(0.1,0.1)
  tau.alpha ~ dgamma(0.1,0.1)
  tau.shift ~ dgamma(0.1,0.1)
  tau.a ~ dgamma(0.1,0.1)
  
}", file="mod_jamil_appli_sensi1.txt")
  
  
  #### Specify the data ####
  sp.data = list(n=n, J=J, S=S, 
                 Z=Z, 
                 alti=alti, survey=c(0,1))
  
  #### Specify the parameters to be monitored ####
  sp.params2 = c('alpha.mean','tau.alpha','alpha',
                 'pmax.mean',
                 'opt.mean','tau.opt', 'opt', 
                 'tol.mean','tau.tol', 'tol', 
                 'shift.mean','tau.shift', 'shift',
                 'tau.a','a.site',
                 'p.fit','p.fitnew'
  )
  
  inits_jamil2 = function() {
    psi.meanGuess=runif(1, 0, 1)
    list(a.site=rnorm(J), 
         pmax.mean=runif(2, 0, 1),
         opt.mean=rnorm(1,1800, 300),
         shift.mean=rnorm(1, 150, 100),
         tol.mean=rnorm(2,330,50))
  }
  
  # MCMC settings
  
  #out_appli_jamil_sensi1 <- autojags(data = sp.data,
  #                                   inits = inits_jamil2,
  #                                   parameters.to.save = sp.params2,
  #                                   model.file = "mod_jamil_appli_sensi1.txt",
  #                                   n.chains = 3,
  #                                   n.adapt = 1000,
  #                                   n.thin = 10,
  #                                   n.burnin = 15000,
  #                                   iter.increment = 15000,
  #                                   DIC = T,
  #                                   store.data = T,
  #                                   Rhat.limit = 1.1,
  #                                   parallel = T,
  #                                   max.iter = 150000)
  
  #save(out_appli_jamil_sensi1, file = "out_appli_jamil_sensi1.RData")
  
  
  #### Model sensitivity 2 : set1 ####
  
  cat("
    model{
  
  for (i in 1:n){
    for (t in 1:S){
      for (j in 1:J){
        Z[j,i,t] ~ dbern(psi.bound[j,i,t])
        psi.bound[j,i,t] <- max(0,min(1,psi[j,i,t]))
        logit(psi[j,i,t]) <- alpha[i,t] - 0.5*pow((alti[j]-(opt[i]+shift[i]*survey[t]))/(tol[i,t]),2) + a.site[j]
      
      #creat simulated dataset to calculate the Bayesain p-value
        Znew[j,i,t] ~ dbern(psi.bound[j,i,t])
        d[j,i,t] <- abs(Z[j,i,t] - psi.bound[j,i,t])
        dnew[j,i,t] <- abs(Znew[j,i,t] - psi.bound[j,i,t])
        d2[j,i,t] <- pow(d[j,i,t],2)
        dnew2[j,i,t] <- pow(dnew[j,i,t],2)
      }
    }
  }
  
  p.fit <- sum(d2[1:J,1:n,1:S])
  p.fitnew <- sum(dnew2[1:J,1:n,1:S])

for (i in 1:n){

    opt[i] ~ dnorm(opt.mean, tau.opt)#T(-4.42, 3.2)
    shift[i] ~ dnorm(shift.mean, tau.shift)#T(-1, 1)

    for (t in 1:S){

    	tol[i,t] ~ dnorm(tol.mean[t], tau.tol)T(0.05,)  
    	alpha[i,t] ~ dnorm(alpha.mean[t], tau.alpha)

    }
  }


for (t in 1:S){
    tol.mean[t] ~ dnorm(330, 1/(70*5)^2) #seems to be way more informative than sensi 1 but in fact the precision parameter (inverse of variance) is still wide for a parameter with the scale of the tolerance                        
    pmax.mean[t] ~ dnorm(0.45,1/(0.27*5)^2)T(0,1)
    alpha.mean[t] <- log(pmax.mean[t])-log(1-pmax.mean[t])
  }
  
    opt.mean ~ dnorm(1800, 1/(240*5)^2) #we can assume that the mean is more likely to be in lower elevations
    shift.mean ~ dnorm(120, 1/(100*5)^2) #we can assume that mean shit is upward, and we could be pretty sure about it

  
  for (j in 1:J){
    a.site[j] ~ dnorm(0, tau.a)
  }
  
  tau.opt ~ dgamma(0.1,0.1)
  tau.tol ~ dgamma(0.1,0.1)
  tau.alpha ~ dgamma(0.1,0.1)
  tau.shift ~ dgamma(0.1,0.1)
  tau.a ~ dgamma(0.1,0.1)
  
}", file="mod_jamil_appli_sensi2.txt")
  
  
  #### Specify the data ####
  sp.data = list(n=n, J=J, S=S, 
                 Z=Z, 
                 alti=alti, survey=c(0,1))
  
  #### Specify the parameters to be monitored ####
  sp.params2 = c('alpha.mean','tau.alpha','alpha',
                 'pmax.mean',
                 'opt.mean','tau.opt', 'opt', 
                 'tol.mean','tau.tol', 'tol', 
                 'shift.mean','tau.shift', 'shift',
                 'tau.a','a.site',
                 'p.fit','p.fitnew'
  )
  
  inits_jamil2 = function() {
    psi.meanGuess=runif(1, 0, 1)
    list(a.site=rnorm(J), 
         pmax.mean=runif(2, 0, 1),
         opt.mean=rnorm(1,1800, 300),
         shift.mean=rnorm(1, 150, 100),
         tol.mean=rnorm(2,330,50))
  }
  
  # MCMC settings
  
  #out_appli_jamil_sensi2 <- autojags(data = sp.data,
  #                                   inits = inits_jamil2,
  #                                   parameters.to.save = sp.params2,
  #                                   model.file = "mod_jamil_appli_sensi2.txt",
  #                                   n.chains = 3,
  #                                   n.adapt = 1000,
  #                                   n.thin = 10,
  #                                   n.burnin = 15000,
  #                                   iter.increment = 15000,
  #                                   DIC = T,
  #                                   store.data = T,
  #                                   Rhat.limit = 1.1,
  #                                   parallel = T,
  #                                   max.iter = 150000)
  
  #save(out_appli_jamil_sensi2, file = "out_appli_jamil_sensi2.RData")
  
  
  
  #### Model sensitivity 2 : set1 ####
  
  cat("
    model{
  
  for (i in 1:n){
    for (t in 1:S){
      for (j in 1:J){
        Z[j,i,t] ~ dbern(psi.bound[j,i,t])
        psi.bound[j,i,t] <- max(0,min(1,psi[j,i,t]))
        logit(psi[j,i,t]) <- alpha[i,t] - 0.5*pow((alti[j]-(opt[i]+shift[i]*survey[t]))/(tol[i,t]),2) + a.site[j]
      
      #creat simulated dataset to calculate the Bayesain p-value
        Znew[j,i,t] ~ dbern(psi.bound[j,i,t])
        d[j,i,t] <- abs(Z[j,i,t] - psi.bound[j,i,t])
        dnew[j,i,t] <- abs(Znew[j,i,t] - psi.bound[j,i,t])
        d2[j,i,t] <- pow(d[j,i,t],2)
        dnew2[j,i,t] <- pow(dnew[j,i,t],2)
      }
    }
  }
  
  p.fit <- sum(d2[1:J,1:n,1:S])
  p.fitnew <- sum(dnew2[1:J,1:n,1:S])

for (i in 1:n){

    opt[i] ~ dnorm(opt.mean, tau.opt)#T(-4.42, 3.2)
    shift[i] ~ dnorm(shift.mean, tau.shift)#T(-1, 1)

    for (t in 1:S){

    	tol[i,t] ~ dnorm(tol.mean[t], tau.tol)T(0.05,)  
    	alpha[i,t] ~ dnorm(alpha.mean[t], tau.alpha)

    }
  }



for (t in 1:S){
    tol.mean[t] ~ dnorm(330, 1/(70*10)^2) #seems to be way more informative than sensi 1 but in fact the precision parameter (inverse of variance) is still wide for a parameter with the scale of the tolerance                        
    pmax.mean[t] ~ dnorm(0.45,1/(0.27*10)^2)T(0,1)
    alpha.mean[t] <- log(pmax.mean[t])-log(1-pmax.mean[t])
  }
  
    opt.mean ~ dnorm(1800, 1/(240*10)^2) #we can assume that the mean is more likely to be in lower elevations
    shift.mean ~ dnorm(120, 1/(100*10)^2) #we can assume that mean shit is upward, and we could be pretty sure about it

  
  for (j in 1:J){
    a.site[j] ~ dnorm(0, tau.a)
  }
  
  tau.opt ~ dgamma(0.1,0.1)
  tau.tol ~ dgamma(0.1,0.1)
  tau.alpha ~ dgamma(0.1,0.1)
  tau.shift ~ dgamma(0.1,0.1)
  tau.a ~ dgamma(0.1,0.1)
  
}", file="mod_jamil_appli_sensi3.txt")
  
  
  #### Specify the data ####
  sp.data = list(n=n, J=J, S=S, 
                 Z=Z, 
                 alti=alti, survey=c(0,1))
  
  #### Specify the parameters to be monitored ####
  sp.params2 = c('alpha.mean','tau.alpha','alpha',
                 'pmax.mean',
                 'opt.mean','tau.opt', 'opt', 
                 'tol.mean','tau.tol', 'tol', 
                 'shift.mean','tau.shift', 'shift',
                 'tau.a','a.site',
                 'p.fit','p.fitnew'
  )
  
  inits_jamil2 = function() {
    psi.meanGuess=runif(1, 0, 1)
    list(a.site=rnorm(J), 
         pmax.mean=runif(2, 0, 1),
         opt.mean=rnorm(1,1800, 300),
         shift.mean=rnorm(1, 150, 100),
         tol.mean=rnorm(2,330,50))
  }
  
  # MCMC settings
  
  #out_appli_jamil_sensi3 <- autojags(data = sp.data,
  #                                   inits = inits_jamil2,
  #                                   parameters.to.save = sp.params2,
  #                                   model.file = "mod_jamil_appli_sensi3.txt",
  #                                   n.chains = 3,
  #                                   n.adapt = 1000,
  #                                   n.thin = 10,
  #                                   n.burnin = 15000,
  #                                   iter.increment = 15000,
  #                                   DIC = T,
  #                                   store.data = T,
  #                                   Rhat.limit = 1.1,
  #                                   parallel = T,
  #                                   max.iter = 150000)
  
  #save(out_appli_jamil_sensi3, file = "out_appli_jamil_sensi3.RData")
}


##############################################%
#### Sensitivity analysis with alti scaled ####
##############################################%

{

#### Model sensitivity 2 : set1 ####

cat("
    model{
  
  for (i in 1:n){
    for (t in 1:S){
      for (j in 1:J){
        Z[j,i,t] ~ dbern(psi.bound[j,i,t])
        psi.bound[j,i,t] <- max(0,min(1,psi[j,i,t]))
        logit(psi[j,i,t]) <- alpha[i,t] - 0.5*pow((alti[j]-(opt[i]+shift[i]*survey[t]))/(tol[i,t]),2) + a.site[j]
      
      #creat simulated dataset to calculate the Bayesain p-value
        Znew[j,i,t] ~ dbern(psi.bound[j,i,t])
        d[j,i,t] <- abs(Z[j,i,t] - psi.bound[j,i,t])
        dnew[j,i,t] <- abs(Znew[j,i,t] - psi.bound[j,i,t])
        d2[j,i,t] <- pow(d[j,i,t],2)
        dnew2[j,i,t] <- pow(dnew[j,i,t],2)
      }
    }
  }
  
  p.fit <- sum(d2[1:J,1:n,1:S])
  p.fitnew <- sum(dnew2[1:J,1:n,1:S])

for (i in 1:n){

    opt[i] ~ dnorm(opt.mean, tau.opt)#T(-4.42, 3.2)
    shift[i] ~ dnorm(shift.mean, tau.shift)#T(-1, 1)

    for (t in 1:S){

    	tol[i,t] ~ dnorm(tol.mean[t], tau.tol)T(0.05,)  
    	alpha[i,t] ~ dnorm(alpha.mean[t], tau.alpha)

    }
  }


for (t in 1:S){
    tol.mean[t] ~ dnorm(0.93, 1/0.20^2) #seems to be way more informative than sensi 1 but in fact the precision parameter (inverse of variance) is still wide for a parameter with the scale of the tolerance                        
    pmax.mean[t] ~ dnorm(0.45,1/0.27^2)T(0,1)
    alpha.mean[t] <- log(pmax.mean[t])-log(1-pmax.mean[t])
  }
  
    opt.mean ~ dnorm(-0.48, 1/0.67^2) #we can assume that the mean is more likely to be in lower elevations
    shift.mean ~ dnorm(0.34, 1/0.28^2) #we can assume that mean shit is upward, and we could be pretty sure about it

  
  for (j in 1:J){
    a.site[j] ~ dnorm(0, tau.a)
  }
  
  tau.opt ~ dgamma(0.1,0.1)
  tau.tol ~ dgamma(0.1,0.1)
  tau.alpha ~ dgamma(0.1,0.1)
  tau.shift ~ dgamma(0.1,0.1)
  tau.a ~ dgamma(0.1,0.1)
  
}", file="mod_jamil_appli_sensi1_scaled.txt")


#### Specify the data ####
sp.data = list(n=n, J=J, S=S, 
               Z=Z, 
               alti=alti.sc, survey=c(0,1))

#### Specify the parameters to be monitored ####
sp.params2 = c('alpha.mean','tau.alpha','alpha',
               'pmax.mean',
               'opt.mean','tau.opt', 'opt', 
               'tol.mean','tau.tol', 'tol', 
               'shift.mean','tau.shift', 'shift',
               'tau.a','a.site',
               'p.fit','p.fitnew'
)

inits_jamil2 = function() {
  psi.meanGuess=runif(1, 0, 1)
  list(a.site=rnorm(J), 
       pmax.mean=runif(2, 0, 1),
       opt.mean=rnorm(1, 0, 0.1),
       shift.mean=rnorm(1, 0, 0.1),
       tol.mean=runif(2,0.1,0.6))
}


# MCMC settings

out_appli_jamil_sensi1_scaled <- autojags(data = sp.data,
                                          inits = inits_jamil2,
                                          parameters.to.save = sp.params2,
                                          model.file = "mod_jamil_appli_sensi1_scaled.txt",
                                          n.chains = 3,
                                          n.adapt = 1000,
                                          n.thin = 10,
                                          n.burnin = 15000,
                                          iter.increment = 15000,
                                          DIC = T,
                                          store.data = T,
                                          Rhat.limit = 1.1,
                                          parallel = T,
                                          max.iter = 150000)

#saveout_appli_jamil_sensi1_scaled, file = "out_appli_jamil_sensi1_scaled.RData")


#### Model sensitivity 2 : set1 ####

cat("
    model{
  
  for (i in 1:n){
    for (t in 1:S){
      for (j in 1:J){
        Z[j,i,t] ~ dbern(psi.bound[j,i,t])
        psi.bound[j,i,t] <- max(0,min(1,psi[j,i,t]))
        logit(psi[j,i,t]) <- alpha[i,t] - 0.5*pow((alti[j]-(opt[i]+shift[i]*survey[t]))/(tol[i,t]),2) + a.site[j]
      
      #creat simulated dataset to calculate the Bayesain p-value
        Znew[j,i,t] ~ dbern(psi.bound[j,i,t])
        d[j,i,t] <- abs(Z[j,i,t] - psi.bound[j,i,t])
        dnew[j,i,t] <- abs(Znew[j,i,t] - psi.bound[j,i,t])
        d2[j,i,t] <- pow(d[j,i,t],2)
        dnew2[j,i,t] <- pow(dnew[j,i,t],2)
      }
    }
  }
  
  p.fit <- sum(d2[1:J,1:n,1:S])
  p.fitnew <- sum(dnew2[1:J,1:n,1:S])

for (i in 1:n){

    opt[i] ~ dnorm(opt.mean, tau.opt)#T(-4.42, 3.2)
    shift[i] ~ dnorm(shift.mean, tau.shift)#T(-1, 1)

    for (t in 1:S){

    	tol[i,t] ~ dnorm(tol.mean[t], tau.tol)T(0.05,)  
    	alpha[i,t] ~ dnorm(alpha.mean[t], tau.alpha)

    }
  }


for (t in 1:S){
    tol.mean[t] ~ dnorm(0.93, 1/(0.20*5)^2) #seems to be way more informative than sensi 1 but in fact the precision parameter (inverse of variance) is still wide for a parameter with the scale of the tolerance                        
    pmax.mean[t] ~ dnorm(0.45,1/(0.27*5)^2)T(0,1)
    alpha.mean[t] <- log(pmax.mean[t])-log(1-pmax.mean[t])
  }
  
    opt.mean ~ dnorm(-0.48, 1/(0.67*5)^2) #we can assume that the mean is more likely to be in lower elevations
    shift.mean ~ dnorm(0.34, 1/(0.28*5)^2) #we can assume that mean shit is upward, and we could be pretty sure about it

  
  for (j in 1:J){
    a.site[j] ~ dnorm(0, tau.a)
  }
  
  tau.opt ~ dgamma(0.1,0.1)
  tau.tol ~ dgamma(0.1,0.1)
  tau.alpha ~ dgamma(0.1,0.1)
  tau.shift ~ dgamma(0.1,0.1)
  tau.a ~ dgamma(0.1,0.1)
  
}", file="mod_jamil_appli_sensi2_scaled.txt")


#### Specify the data ####
sp.data = list(n=n, J=J, S=S, 
               Z=Z, 
               alti=alti.sc, survey=c(0,1))

#### Specify the parameters to be monitored ####
sp.params2 = c('alpha.mean','tau.alpha','alpha',
               'pmax.mean',
               'opt.mean','tau.opt', 'opt', 
               'tol.mean','tau.tol', 'tol', 
               'shift.mean','tau.shift', 'shift',
               'tau.a','a.site',
               'p.fit','p.fitnew'
)

inits_jamil2 = function() {
  psi.meanGuess=runif(1, 0, 1)
  list(a.site=rnorm(J), 
       pmax.mean=runif(2, 0, 1),
       opt.mean=rnorm(1, 0, 0.1),
       shift.mean=rnorm(1, 0, 0.1),
       tol.mean=runif(2,0.1,0.6))
}


# MCMC settings

out_appli_jamil_sensi2_scaled <- autojags(data = sp.data,
                                          inits = inits_jamil2,
                                          parameters.to.save = sp.params2,
                                          model.file = "mod_jamil_appli_sensi2_scaled.txt",
                                          n.chains = 3,
                                          n.adapt = 1000,
                                          n.thin = 10,
                                          n.burnin = 15000,
                                          iter.increment = 15000,
                                          DIC = T,
                                          store.data = T,
                                          Rhat.limit = 1.1,
                                          parallel = T,
                                          max.iter = 150000)

#saveout_appli_jamil_sensi2_scaled, file = "out_appli_jamil_sensi2_scaled.RData")



#### Model sensitivity 2 : set1 ####

cat("
    model{
  
  for (i in 1:n){
    for (t in 1:S){
      for (j in 1:J){
        Z[j,i,t] ~ dbern(psi.bound[j,i,t])
        psi.bound[j,i,t] <- max(0,min(1,psi[j,i,t]))
        logit(psi[j,i,t]) <- alpha[i,t] - 0.5*pow((alti[j]-(opt[i]+shift[i]*survey[t]))/(tol[i,t]),2) + a.site[j]
      
      #creat simulated dataset to calculate the Bayesain p-value
        Znew[j,i,t] ~ dbern(psi.bound[j,i,t])
        d[j,i,t] <- abs(Z[j,i,t] - psi.bound[j,i,t])
        dnew[j,i,t] <- abs(Znew[j,i,t] - psi.bound[j,i,t])
        d2[j,i,t] <- pow(d[j,i,t],2)
        dnew2[j,i,t] <- pow(dnew[j,i,t],2)
      }
    }
  }
  
  p.fit <- sum(d2[1:J,1:n,1:S])
  p.fitnew <- sum(dnew2[1:J,1:n,1:S])

for (i in 1:n){

    opt[i] ~ dnorm(opt.mean, tau.opt)#T(-4.42, 3.2)
    shift[i] ~ dnorm(shift.mean, tau.shift)#T(-1, 1)

    for (t in 1:S){

    	tol[i,t] ~ dnorm(tol.mean[t], tau.tol)T(0.05,)  
    	alpha[i,t] ~ dnorm(alpha.mean[t], tau.alpha)

    }
  }



for (t in 1:S){
    tol.mean[t] ~ dnorm(0.93, 1/(0.20*10)^2) #seems to be way more informative than sensi 1 but in fact the precision parameter (inverse of variance) is still wide for a parameter with the scale of the tolerance                        
    pmax.mean[t] ~ dnorm(0.45,1/(0.27*10)^2)T(0,1)
    alpha.mean[t] <- log(pmax.mean[t])-log(1-pmax.mean[t])
  }
  
    opt.mean ~ dnorm(-0.48, 1/(0.67*10)^2) #we can assume that the mean is more likely to be in lower elevations
    shift.mean ~ dnorm(0.34, 1/(0.28*10)^2) #we can assume that mean shit is upward, and we could be pretty sure about it

  
  for (j in 1:J){
    a.site[j] ~ dnorm(0, tau.a)
  }
  
  tau.opt ~ dgamma(0.1,0.1)
  tau.tol ~ dgamma(0.1,0.1)
  tau.alpha ~ dgamma(0.1,0.1)
  tau.shift ~ dgamma(0.1,0.1)
  tau.a ~ dgamma(0.1,0.1)
  
}", file="mod_jamil_appli_sensi3_scaled.txt")


#### Specify the data ####
sp.data = list(n=n, J=J, S=S, 
               Z=Z, 
               alti=alti.sc, survey=c(0,1))

#### Specify the parameters to be monitored ####
sp.params2 = c('alpha.mean','tau.alpha','alpha',
               'pmax.mean',
               'opt.mean','tau.opt', 'opt', 
               'tol.mean','tau.tol', 'tol', 
               'shift.mean','tau.shift', 'shift',
               'tau.a','a.site',
               'p.fit','p.fitnew'
)

inits_jamil2 = function() {
  psi.meanGuess=runif(1, 0, 1)
  list(a.site=rnorm(J), 
       pmax.mean=runif(2, 0, 1),
       opt.mean=rnorm(1, 0, 10),
       shift.mean=rnorm(1, 0, 10),
       tol.mean=rnorm(2, 50, 10))
}

# MCMC settings

out_appli_jamil_sensi3_scaled <- autojags(data = sp.data,
                                          inits = inits_jamil2,
                                          parameters.to.save = sp.params2,
                                          model.file = "mod_jamil_appli_sensi3_scaled.txt",
                                          n.chains = 3,
                                          n.adapt = 1000,
                                          n.thin = 10,
                                          n.burnin = 15000,
                                          iter.increment = 15000,
                                          DIC = T,
                                          store.data = T,
                                          Rhat.limit = 1.1,
                                          parallel = T,
                                          max.iter = 150000)

#save(out_appli_jamil_sensi3_scaled, file = "out_appli_jamil_sensi3_scaled.RData")  
}



save(out_appli_jamil_sensi1,
     out_appli_jamil_sensi2,
     out_appli_jamil_sensi3,
     out_appli_jamil_sensi1_scaled,
     out_appli_jamil_sensi2_scaled,
     out_appli_jamil_sensi3_scaled,
     file = "~/Bureau/PhD/Orthoptera/CC/Analysis/Application_analysis/Results_data/out_appli_sensitivity_analysis.RData")
