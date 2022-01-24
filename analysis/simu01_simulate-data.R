############################################%
# Simulate data for EHMOS paper
#
# last modification: 09/12/21 (edit annotations)
# bastien.mourguiart@gmail.com
#
# SIMU_data-simulation.R
#
# Script to create virtual presence/absence data
# used in the simulation study in the paper about optimum shifts.
#
# Simulate eight scenarios that result of the combination of three sub scenarios
# each decomposed in two categories.
# Outputs are true ecological parameters for each species (especially true shifts),
# one occupancy data set per scenario, and 30 pres/abs dataset per scenario
############################################%

set.seed(1904) #set a seed for replicability

## Create data frame that will contain simulated parameters ##
## Set the number of sites and species ##
J <- 300
N <- 20
sim.data <- data.frame("scenario" = rep(c("A1xB1xE1", # scenario names had evolved during the work, E1 became C1 after
                                          "A1xB2xE1", # at the start E1 refers to Ecological sub-scenario 1
                                          "A1xB1xE2", # but became C1 for simplicity and match with A & B
                                          "A1xB2xE2",
                                          "A2xB1xE1",
                                          "A2xB2xE1",
                                          "A2xB1xE2",
                                          "A2xB2xE2"
                                          ),
                                        each = N
                                        ),
                       "sp" = paste("sp",1:N, sep = ""),
                       "spType" = NA,
                       "T_opti" = NA,
                       "T_shift" = NA,
                       "T_shape" = NA,
                       "Opti" = NA,
                       "Shift" = NA,
                       "Width1" = NA,
                       "Pmax1" = NA,
                       "Width2" = NA,
                       "Pmax2" = NA,
                       "Tol1" = NA,
                       "Alpha1" = NA,
                       "Tol2" = NA,
                       "Alpha2" = NA)

#### 1. Simulate sampling sites ####

## Distribute the sites on the gradient according to the scenario ##
#Scenario with normal sampling
X_A1 <- msm::rmenorm(J, 1970, 355, 1000, 3000) #mean and sd chose accordingly with the application design
#Scenario with uniform sampling
X_A2 <- runif(J, 1000, 3000)

## Deduce the gradient 'limits' in which optimum could be simulated depending on sampling design ##
L_minA1 <- quantile(X_A1, 0.02) # below the bottom edge, 10% of sites remain
L_lowA1 <- quantile(X_A1, 0.1) # set a minimum for optimum within the sampling range
L_highA1 <- quantile(X_A1, 0.9)
L_maxA1 <- quantile(X_A1, 0.98)
L_minA2 <- quantile(X_A2, 0.02)
L_lowA2 <- quantile(X_A2, 0.1)
L_highA2 <- quantile(X_A2, 0.9)
L_maxA2 <- quantile(X_A2, 0.98)

#### 2. Simulate ecological parameters ####
#### 2.1. Simulate optima and shifts ####

## Create matrix containing species type depending on sub-scenario combinations
mat.NBtype <- matrix(data = c(0, 0, 0, 0, 20, 0,  #B1xE1,  only middle and specialist species
                              0, 0, 3, 3, 14, 0,  #B2xE1,  14 middle and 6 edge (3 at each boundary),  all specialist
                              0, 0, 0, 0, 10, 10,  #B1xE2,  all middle,  50-50 specialist/generalist
                              2, 1, 1, 2, 7, 7 #B2xE2,  14 middle and 6 edge, 50-50 specialist/generalist
                              ),
                     nrow = 4,
                     byrow = T,
                     dimnames = list(c("B1xE1","B2xE1","B1xE2","B2xE2"),
                                     c("Eup_G","Edn_G","Eup_S","Edn_S","Mid_S","Mid_G")))

library(foreach)
foreach::foreach (s = unique(sim.data$scenario)) %do% {
  if(substr(s, 1, 2)=="A1"){
    X <- X_A1

    ## Deduce the sampling 'limits' for each scenario ##
    L_min <- L_minA1
    L_low <- L_lowA1
    L_high <- L_highA1
    L_max <- L_maxA1
  }else{
    X <- X_A2

    ## Deduce the sampling 'limits' for each scenario ##
    L_min <- L_minA2
    L_low <- L_lowA2
    L_high <- L_highA2
    L_max <- L_maxA2
  }

  if(substr(s, 4, 8) == "B1xE1"){ # set the expected proportion of edge species depending on the subscenario
    NBtype <- as.numeric(mat.NBtype[1,])
    e <- 0 # proportion of edge species in subscenario b1
  }
  if(substr(s, 4, 8) == "B2xE1"){
    NBtype <- as.numeric(mat.NBtype[2,])
    e <- 0.3 # proportion of edge species in subscenario b1
  }
  if(substr(s, 4, 8) == "B1xE2"){
    NBtype <- as.numeric(mat.NBtype[3,])
    e <- 0
  }
  if(substr(s, 4, 8) == "B2xE2"){
    NBtype <- as.numeric(mat.NBtype[4,])
    e <- 0.3
  }

  sim.data[sim.data$scenario == s,]$spType <- sample(x = rep(colnames(mat.NBtype), NBtype), size = 20, replace = F) # attribute randomly the species type depending on the scenario
  sim.data[sim.data$scenario == s,]$T_opti <- substr(sim.data[sim.data$scenario == s,]$spType, 1, 3) # deduce the species optimum type
  sim.data[sim.data$scenario == s,]$T_shift <- sample(x=rep(c("Up","Down","None"), c(12,4,4)), size = 20, replace = F) # attribute randomly the type of shift, with same proportion of types for all scenarios
  sim.data[sim.data$scenario == s,]$T_shape <- substr(sim.data[sim.data$scenario == s,]$spType, 5, 5) # deduce the species specialization type

  h<-1 ; while(h<2){
    for (i in 1:N){ #for each species set an optimum and a shift depending on its types
      if (sim.data[sim.data$scenario == s,]$T_opti[i] == "Eup"){
        sim.data[sim.data$scenario == s,]$Opti[i] <- runif(1,L_high,L_max)
      }
      if (sim.data[sim.data$scenario == s,]$T_opti[i] == "Edn"){
        sim.data[sim.data$scenario == s,]$Opti[i] <- runif(1,L_min,L_low)
      }
      if (sim.data[sim.data$scenario == s,]$T_opti[i] == "Mid"){
        sim.data[sim.data$scenario == s,]$Opti[i] <- runif(1,L_low+200,L_high-200)
      }
      if (sim.data[sim.data$scenario == s,]$T_shift[i] == "Up"){
        sim.data[sim.data$scenario == s,]$Shift[i] <- runif(1, 80, 200)
      }
      if (sim.data[sim.data$scenario == s,]$T_shift[i] == "Down"){
        sim.data[sim.data$scenario == s,]$Shift[i] <- runif(1, -200, -80)
      }
      if (sim.data[sim.data$scenario == s,]$T_shift[i] == "None"){
        sim.data[sim.data$scenario == s,]$Shift[i] <- 0
      }
    }

    #Verify if the number of edge species is good
    # Needed because depending on shift simulation, middle species could become edge species in second sampling occasion which is not wanted
    Edg<- ((sim.data[sim.data$scenario == s,]$Opti > L_high & (sim.data[sim.data$scenario == s,]$Opti + sim.data[sim.data$scenario == s,]$Shift) > L_high) | (sim.data[sim.data$scenario == s,]$Opti < L_low & (sim.data[sim.data$scenario == s,]$Opti + sim.data[sim.data$scenario == s,]$Shift) < L_low))
    if(sum(Edg)==(e*N)){ # verify if the simulated number of edge species corresponds to what expected
      h <- 2
    }
  }

  #### 2.2. Simulate max probabilities and widths ####
  Pmax <- matrix(0, nrow=2, ncol=N)
  w <- matrix(0, nrow=2, ncol=N)

  for (i in 1:N){ # simulate max probability and width according to subscenario
    if (sim.data[sim.data$scenario == s,]$T_shape[i]=="S"){
      sim.data[sim.data$scenario == s,]$Pmax1[i] <- runif(1, 0.9, 0.99)
      sim.data[sim.data$scenario == s,]$Pmax2[i] <- runif(1, 0.9, 0.99) # values could vary between sampling occasions
      sim.data[sim.data$scenario == s,]$Width1[i] <- runif(1, 400, 600)
      sim.data[sim.data$scenario == s,]$Width2[i] <- runif(1, 400, 600)
    }else{
      sim.data[sim.data$scenario == s,]$Pmax1[i] <- runif(1, 0.75, 0.85)
      sim.data[sim.data$scenario == s,]$Pmax2[i] <- runif(1, 0.75, 0.85)
      sim.data[sim.data$scenario == s,]$Width1[i] <- runif(1, 1200, 1400)
      sim.data[sim.data$scenario == s,]$Width2[i] <- runif(1, 1200, 1400)
    }
  }
}

#### 3. Transform ecological parameters into model parameters ####
sim.data$Alpha1 <- log(sim.data$Pmax1/(1-sim.data$Pmax1))
sim.data$Alpha2 <- log(sim.data$Pmax2/(1-sim.data$Pmax2))
sim.data$Tol1 <- 0.5*sim.data$Width1/sqrt(2*(sim.data$Alpha1 - log(0.05/0.95)))
sim.data$Tol2 <- 0.5*sim.data$Width2/sqrt(2*(sim.data$Alpha2 - log(0.05/0.95)))

#### 4. Deduce from parameters species-environment relationship ####
## Create one unique array containing simulated occupancy probabilities for each species in each scenarios ##
psi <- array(0, dim=c(J,N,2,8),
             dimnames = list(paste("site",1:J, sep=""),
                             paste("sp",1:N, sep=""),
                             c("p1","p2"),
                             c("A1xB1xE1","A1xB2xE1","A1xB1xE2","A1xB2xE2","A2xB1xE1","A2xB2xE1","A2xB1xE2","A2xB2xE2")))

for (s in c("A1xB1xE1","A1xB2xE1","A1xB1xE2","A1xB2xE2","A2xB1xE1","A2xB2xE1","A2xB1xE2","A2xB2xE2")){
  if(substr(s, 1, 2)=="A1"){
    X <- X_A1
  }else{
    X <- X_A2
  }
  for (i in 1:N){
    psi[,i,1,s] <- plogis(sim.data[sim.data$scenario == s,]$Alpha1[i] - ((X - sim.data[sim.data$scenario == s,]$Opti[i])^2/(2*sim.data[sim.data$scenario == s,]$Tol1[i]^2)))
    psi[,i,2,s] <- plogis(sim.data[sim.data$scenario == s,]$Alpha2[i] - ((X - (sim.data[sim.data$scenario == s,]$Opti[i]+sim.data[sim.data$scenario == s,]$Shift[i]))^2/(2*sim.data[sim.data$scenario == s,]$Tol2[i]^2)))
  }
}

#### 5. Produce replicated presence absence datasets from species-environment relationships ####
## Create 30 replicated data sets of presence absence ##
Z <- array(0, dim=c(J,N,2,8,30),
           dimnames = list(paste("site",1:J, sep=""),
                           paste("sp",1:N, sep=""),
                           c("p1","p2"),
                           c("A1xB1xE1","A1xB2xE1","A1xB1xE2","A1xB2xE2","A2xB1xE1","A2xB2xE1","A2xB1xE2","A2xB2xE2"),
                           paste("repli",1:30, sep="")))

for (s in c("A1xB1xE1","A1xB2xE1","A1xB1xE2","A1xB2xE2","A2xB1xE1","A2xB2xE1","A2xB1xE2","A2xB2xE2")){
  for (p in 1:2){
    for (i in 1:N){
      for (j in 1:J){
        Z[j,i,p,s,] <- rbinom(30, 1, psi[j,i,p,s])
      }
    }
  }
}

# Simplify opti types to 2 opti type (mid, edge) in sim.data
sim.data %>%
  dplyr::mutate("T_opti" = forcats::as_factor(T_opti)
  ) %>%
  dplyr::mutate("T_opti" = forcats::fct_collapse(T_opti,
                                                 Edge = c("Eup", "Edn"))
                ) -> sim.data

save(list = c("sim.data", "psi", "X_A1", "X_A2"),
     file = here::here("results", "simulated_data.RData"))


