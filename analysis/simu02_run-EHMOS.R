# R code for the simulation study in:
# A new method to explicitly estimate the shift of optimum along gradients in multispecies studies.
# B. Mourguiart, B. Liquet, K. Mengersen, T. Couturier, J. Mansons, Y. Braud, A. Besnard

# Script 2: Fit EHMOS model to simulated data

# This script can take several days to run
# Models' outputs created by this script ("simu_out-ehmos.RData") 
# can be download at: https://doi.org/10.5061/dryad.0gb5mkm59

load(here::here("results", "simulated_data.RData")) # created in simu01_simulate-data.R

## Scale X
m.site_A1 <- mean(X_A1)
sd.site_A1 <- sd(X_A1)
sc.site_A1 <- (X_A1 - m.site_A1) / sd.site_A1

m.site_A2 <- mean(X_A2)
sd.site_A2 <- sd(X_A2)
sc.site_A2 <- (X_A2 - m.site_A2) / sd.site_A2


#### cat model EHMOS multi sp ####
#use model with random alpha and tolerance

cat("
    model{

  for (i in 1:n){
    for (t in 1:S){
      for (j in 1:J){
        z[j,i,t] ~ dbern(psi.bound[j,i,t])
        psi.bound[j,i,t] <- max(0,min(1,psi[j,i,t]))
        logit(psi[j,i,t]) <- b0[i,t] - 0.5*pow((alti[j]-(opt[i]+shift[i]*survey[t]))/(tol[i,t]),2) + a.site[j]

      #creat simulated dataset to calculate the Bayesain p-value
        znew[j,i,t] ~ dbern(psi.bound[j,i,t])
        d[j,i,t] <- abs(z[j,i,t] - psi.bound[j,i,t])
        dnew[j,i,t] <- abs(znew[j,i,t] - psi.bound[j,i,t])
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
    	b0[i,t] ~ dnorm(b0.mean[t], tau.b0)

    }
  }


for (t in 1:S){
    tol.mean[t] ~ dnorm(0, 0.001)T(0.05,)
    b0.mean[t] ~ dnorm(0, 0.001)
  }

    opt.mean ~ dnorm(0, 0.001)
    shift.mean ~ dnorm(0, 0.001)


  for (j in 1:J){
    a.site[j] ~ dnorm(0, tau.a)
  }

  tau.opt ~ dgamma(0.1,0.1)
  tau.tol ~ dgamma(0.1,0.1)
  tau.b0 ~ dgamma(0.1,0.1)
  tau.shift ~ dgamma(0.1,0.1)
  tau.a ~ dgamma(0.1,0.1)

}", file="mod_EHMOS.txt")

inits_EHMOS = function() {
  psi.meanGuess=runif(1, 0, 1)
  list(a.site=rnorm(300),
       b0.mean=rnorm(2, mean=2, sd=1),
       opt.mean=runif(1,-1.6, 0.1),
       shift.mean=rnorm(1, 0, 0.1),
       tol.mean=runif(2,0.1,0.6))
}

#### Specify the parameters to be monitored
params_EHMOS = c('opt.mean','tau.opt', 'opt',
                 'tol.mean','tau.tol', 'tol',
                 'b0.mean','tau.b0', 'b0',
                 'shift.mean','tau.shift', 'shift',
                 'tau.a','a.site',
                 'p.fit','p.fitnew'
)

library(doParallel)
cl <- doParallel::makeForkCluster(30)
doParallel::registerDoParallel(cl)


res.ehmos <- foreach (s=1:8) %:%
	foreach (r=1:30) %dopar% {
	    set.seed(64)
	      data_jags <- list(n=20, J=300, S=2,
		                z=Z[,,,s,r],
		                alti=sc.site_A2,
		                survey=c(0,1))
	      name <- paste("out.sim2_", dimnames(Z)[[4]][s],"_EHMOS", r, sep="")
	      out <- autojags(data = data_jags,
		          inits = inits_EHMOS,
		          parameters.to.save = params_EHMOS,
		          model.file = "mod_EHMOS.txt",
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
	      out
  }
save(res.ehmos, file=here::here("results", "simu_out-ehmos.RData"))

