############################################%
# Fit Bayesian GLMM on simulated data
#
# last modification: 09/12/21 (edit annotations)
# bastien.mourguiart@gmail.com
#
# simu02_run-GLMM.R
# depends on: simu01_simulate-data.R
#
# Script coding and running a Bayesian GLMM
# on simulated data.
# Script to run in a cluster, or segmented.
# Take several days to run.
############################################%


load(here::here("results", "simulated_data.RData")) # created in simu01_simulate-data.R


## Scale X
m.site_A1 <- mean(X_A1)
sd.site_A1 <- sd(X_A1)
sc.site_A1 <- (X_A1 - m.site_A1) / sd.site_A1

m.site_A2 <- mean(X_A2)
sd.site_A2 <- sd(X_A2)
sc.site_A2 <- (X_A2 - m.site_A2) / sd.site_A2

cat("
    model{

    #Define prior distributions for community-level model parameters
    #b0.mean ~ dunif(0,1)                          #species level effect on occupancy probability
    #mu.b0 <- log(b0.mean) - log(1-b0.mean)        #species level effect on occupancy probability on the logit scale

    mu.b0 ~ dnorm(0, 0.001)
    mub1 ~ dnorm(0, 0.001)                        #detection covariable 1 effect on detectability
    mub2 ~ dnorm(0, 0.001)
    mub3 ~ dnorm(0, 0.001)
    mub4 ~ dnorm(0, 0.001)
    mub5 ~ dnorm(0, 0.001)

    tau.a ~ dgamma(0.1,0.1)                      #variability of species level effect on detectability
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
            z[j,i,t] ~ dbern(psi.bound[j,i,t])
            psi.bound[j,i,t]<-max(0,min(1,psi[j,i,t]))
            logit(psi[j,i,t]) <- b0[i] + b1[i]*alti[j] + b2[i]*pow(alti[j],2) +
                                 b3[i]*survey[t] + b4[i]*alti[j]*survey[t] + b5[i]*pow(alti[j],2)*survey[t] +
        		                     a[j]
            #creat simulated dataset to calculate the Bayesain p-value
        Znew[j,i,t] ~ dbern(psi.bound[j,i,t])
        d[j,i,t] <- abs(z[j,i,t] - psi.bound[j,i,t])
        dnew[j,i,t] <- abs(Znew[j,i,t] - psi.bound[j,i,t])
        d2[j,i,t] <- pow(d[j,i,t],2)
        dnew2[j,i,t] <- pow(dnew[j,i,t],2)
      }
    }
  }

  p.fit <- sum(d2[1:J,1:n,1:S])
  p.fitnew <- sum(dnew2[1:J,1:n,1:S])
    }

", file="mod_glmm.txt")


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
  list(a=rnorm(300),
       mu.b0=rnorm(1),
       mub1=rnorm(1),
       mub2=rnorm(1),
       mub3=rnorm(1),
       mub4=rnorm(1),
       mub5=rnorm(1))
}



library(doParallel)
cl <- makeForkCluster(30)
registerDoParallel(cl)


res.glmm <- foreach::foreach (s=1:8) %:%
  foreach::foreach (r=1:30) %dopar% {
    set.seed(64) # maybe useless
    data_jags <- list(n=20, J=300, S=2,
                      z=Z[,,,s,r],
                      alti=sc.site_A2,
                      survey=c(0,1))
    name <- paste("out.sim2_", dimnames(Z)[[4]][s],"_glmm", r, sep="")
    out <- jagsUI::autojags(data = data_jags,
                            inits = inits_glmm,
                            parameters.to.save = params_glmm,
                            model.file = "mod_glmm.txt",
                            n.chains = 3,
                            n.adapt = 500,
                            n.thin = 15,
                            n.burnin = 15000,
                            iter.increment = 15000,
                            DIC = TRUE,
                            store.data = TRUE, #keep initial values
                            parallel = F,
                            max.iter = 250000)
    save(out, file = here::here("results", paste0(name, ".RData"))) # in case of crash
  }

save(res.glmm, file=here::here("results", "simu_out-glmm.RData"))



