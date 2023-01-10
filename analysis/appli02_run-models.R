# R code for the application study in:
# A new method to explicitly estimate the shift of optimum along gradients in multispecies studies.
# B. Mourguiart, B. Liquet, K. Mengersen, T. Couturier, J. Mansons, Y. Braud, A. Besnard

# Script 2: Fit Bayesian models to Orthoptera data

set.seed(44)

alti.sc <- as.numeric(scale(alti, center = TRUE, scale = TRUE))

#### 1. EHMOS ----

## .. 1.1. Write the model ####
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
    p.fitSP[i] <- sum(d2[1:J,i,1:S])
    p.fitnewSP[i] <- sum(dnew2[1:J,i,1:S])
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

}", file="mod_EHMOS.txt")

## .. 1.2 Specify the data ####
sp.data <- list(n = n,
                J = J,
                S = S,
                Z = Z,
                alti = alti.sc,
                survey=c(0,1)
                )

## .. 1.3 Specify the parameters to be monitored ####
sp.params <- c('alpha.mean', 'tau.alpha', 'alpha',
                'opt.mean', 'tau.opt',  'opt',
                'tol.mean', 'tau.tol',  'tol',
                'shift.mean', 'tau.shift',  'shift',
                'tau.a', 'a.site',
                'p.fit', 'p.fitnew',
                'p.fitSP', 'p.fitnewSP'
               )

## .. 1.4. Specify intial values ####
inits_ehmos <- function() {
  psi.meanGuess = runif(1, 0, 1)
  list(a.site = rnorm(J),
       alpha.mean = rnorm(2, mean = 2, sd = 1),
       opt.mean = rnorm(1,1800, 300),
       shift.mean = rnorm(1, 150, 100),
       tol.mean = rnorm(2, 330, 50))
}

## .. 1.5. Run the model ####
out_appli_ehmos <- jagsUI::autojags(data = sp.data,
                                    inits = inits_ehmos,
                                    parameters.to.save = sp.params,
                                    model.file = "mod_EHMOS.txt",
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
## .. 1.6. Save outputs ####
save(out_appli_ehmos,
     file = here::here("results", "appli_out-ehmos.RData"))



#### 2. cGLMM ----
## 2.1. Write the model ####
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
    p.fitSP[i] <- sum(d2[1:J,i,1:S])
    p.fitnewSP[i] <- sum(dnew2[1:J,i,1:S])
  }

  p.fit <- sum(d2[1:J,1:n,1:S])
  p.fitnew <- sum(dnew2[1:J,1:n,1:S])
    }

", file = "mod_glmm.txt")

## .. 2.2. Specify parameters to save ####
params_glmm <- c('mu.b0', 'mub1', 'mub2', 'mub3', 'mub4', 'mub5',
                 'tau.b0', 'tau.b1', 'tau.b2', 'tau.b3', 'tau.b4', 'tau.b5',
                 'tau.a',
                 'b0', 'b1', 'b2', 'b3', 'b4', 'b5',
                 'a',
                 'p.fit', 'p.fitnew',
                 'p.fitSP',  'p.fitnewSP'
)

## .. 2.3. Specify initial values ####
inits_glmm <- function() {
  psi.meanGuess = runif(1, 0, 1)
  list(a = rnorm(J),
       b0.mean = runif(1, 0, 1),
       mub1 = rnorm(1),
       mub2 = rnorm(1),
       mub3 = rnorm(1),
       mub4 = rnorm(1),
       mub5 = rnorm(1))
}

## .. 2.4. Run the model ####
out_appli_glmm <- jagsUI::autojags(data = sp.data,
                                   inits = inits_glmm,
                                   parameters.to.save = params_glmm,
                                   model.file = "mod_glmm.txt",
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

## .. 2.5. Save outputs ####
save(out_appli_glmm,
     file = here::here("results", "appli_out-glmm.RData"))
