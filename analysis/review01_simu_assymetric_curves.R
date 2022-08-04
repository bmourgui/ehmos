
library(magrittr)
library(ggplot2)

set.seed(1904)


#### Sampling scenario: same as first simulations ####
J <- 300 # number of sampling sites
N <- 9 # number of species
X_A1 <- msm::rmenorm(J, 1970, 355, 1000, 3000) #mean and sd chose accordingly with the application design
X_A1 <- X_A1[order(X_A1)]
#Scenario with uniform sampling
X_A2 <- runif(J, 1000, 3000)
X_A2 <- X_A2[order(X_A2)]

## Deduce the gradient 'limits' in which optimum could be simulated depending on sampling design ##
L_minA1 <- quantile(X_A1, 0.02) # below the bottom edge, 10% of sites remain
L_lowA1 <- quantile(X_A1, 0.1) # set a minimum for optimum within the sampling range
L_highA1 <- quantile(X_A1, 0.9)
L_maxA1 <- quantile(X_A1, 0.98)
L_minA2 <- quantile(X_A2, 0.02)
L_lowA2 <- quantile(X_A2, 0.1)
L_highA2 <- quantile(X_A2, 0.9)
L_maxA2 <- quantile(X_A2, 0.98)

#### Species Response Curves: new ####

# Create a function to simulate species response curves.
# Simulated SRCs are mixtures of two normal distributions having the same mean,
# but with potentially different variances. First distribution (resp. second) is used to 
# simulate the left (resp. right) side of optimum (i.e. the mean of the distribution) of the SRC.
# Hence, if the first distribution has a greater variance than the second distribution,
# the curve will be asymmetric to the left. If both distribution have the same variance, the
# SRC is symmetric. Note that we used the density function of the normal distribution, but
# with a parameter scale to reach a specified maximum probability value.
src_asym <- function(x, opt, sd_l, sd_r, max){
  y <- exp(-(x - opt)^2 / (2 * (sd_l ^ 2))) / (sd_l * sqrt(2 * 3.14))
  y2 <- exp(-(x - opt)^2 / (2 * (sd_r ^ 2))) / (sd_r * sqrt(2 * 3.14))
  y <- y * (max / max(y))
  y2 <- y2 * (max / max(y2))
  
  xy <- data.frame(x = x, y = y, y2 = y2)
  
  c(xy[xy$x <= opt, ]$y, xy[xy$x > opt, ]$y2)
}

## determine widths of the curve at psi = 0.05 for different sd 
qnorm(0.95, mean = 1900, sd = seq(110, 200, 10)) - qnorm(0.05, mean = 1900, sd = seq(110, 200, 10))
# -> if symmetric widths between 400 and 600 m (specialist width) correspond to (approx.) sd between 125 and 185
(qnorm(0.95, mean = 1900, sd = 120) - qnorm(0.5, mean = 1900, sd = 120)) +  (qnorm(0.5, mean = 1900, sd = 170) - qnorm(0.05, mean = 1900, sd = 170))
# -> asymmetric one curve with sd 120 and the other with 170 has a width of 477 m
sd2 <- rnorm(1000, 130, 10)
sd1 <- rnorm(1000, 180, 10)
range((qnorm(0.95, mean = 1900, sd = sd2) - qnorm(0.5, mean = 1900, sd = sd2)) +  (qnorm(0.5, mean = 1900, sd = sd1) - qnorm(0.05, mean = 1900, sd = sd1)))
# -> if we allow variation between species


# We defined three type of potential species response curves: symmetric, asymmetric to the left
# and to the right. We simulated 9 species to describe all possible combinations of SRC that
# could be between two periods. To keep it simple we chose to only simulate specialist species
# with optimum positions in the middle of the sampling range according to sampling scenario A1.
# Besides, we randomly attribute upward shifts following random draws from a normal distribution
# with mean and standard deviation of  100 m and 10 m respectively. 
# Once SRCs were simulated for each species and periods, we virtually observed
# presence/absence data from random draws of a Bernoulli distribution with parameter gave by SRCs. 
# We replicated the process 30 times.

# create data with species attributes used to calculate SRC
data.sp <- data.frame(
  # define type of curves for period 1
  "curve.type1" = rep(c("sym", #symmetric
                        "asym_r", #asymmetric to the right
                        "asym_l"), #asymmetric to the left
                      3),
  # define type of curves (toc) for period 1
  "curve.type2" = rep(c("sym", #symmetric
                        "asym_r", #asymmetric to the right
                        "asym_l"), #asymmetric to the left
                      each = 3),
  # according to the toc, define sd of the right and lefts parts of normal curves
  "sd_l.1" = rep(c(100, 100, 150), 3), # sd left part first period
  "sd_r.1" = rep(c(100, 150, 100), 3), # sd right part first period
  "sd_l.2" = rep(c(100, 100, 150), each = 3), # sd left part second period
  "sd_r.2" = rep(c(100, 150, 100), each = 3),# sd right part second period
  # simulate optima with only middle species
  "opt" = runif(9, L_lowA1 + 200, L_highA1 - 200),
  # simulate only upward shifts
  "shift" = runif(9, 80, 120),
  # simulate only specialist species
  "max.occ1" = runif(9, 0.85, 0.99),
  "max.occ2" = runif(9, 0.85, 0.99)
)


## Deduce from parameters species-environment relationship 
## Create one unique array containing simulated occupancy probabilities for each species in each scenarios ##
psi <- array(0, dim=c(J, N, 2, 2),
             dimnames = list(paste("site", 1:J, sep=""),
                             paste("sp", 1:N, sep=""),
                             c("p1", "p2"),
                             c("A1", "A2")))

for (s in c("A1", "A2")){
  if(s == "A1"){
    X <- X_A1
  }else{
    X <- X_A2
  }
  for (i in 1:N){
    psi[,i,1,s] <- src_asym(x = X[order(X)], 
                            opt = data.sp[i, ]$opt, 
                            sd_l = data.sp[i, ]$sd_l.1, 
                            sd_r = data.sp[i, ]$sd_r.1,
                            max = data.sp[i, ]$max.occ1)
    psi[,i,2,s] <- src_asym(x = X[order(X)], 
                            opt = (data.sp[i, ]$opt + data.sp[i, ]$shift), 
                            sd_l = data.sp[i, ]$sd_l.2, 
                            sd_r = data.sp[i, ]$sd_r.2,
                            max = data.sp[i, ]$max.occ2)
  }
}

## Create 30 replicated data sets of presence absence ##
Z <- array(0, dim=c(J, N, 2, 2, 30),
           dimnames = list(paste("site",1:J, sep=""),
                           paste("sp",1:N, sep=""),
                           c("p1","p2"),
                           c("A1", "A2"),
                           paste("repli", 1:30, sep="")))

for (s in c("A1", "A2")){
  for (p in 1:2){
    for (i in 1:N){
      for (j in 1:J){
        Z[j, i, p, s, ] <- rbinom(30, 1, psi[j, i, p, s])
      }
    }
  }
}

#### Visualization ####

# Illustration with three species for one period:
# symmetric (black), asymmetric left (red) and right (green)
x <- 1000:3000
y.sym <- src_asym(x = x, opt = 2000, sd_l = 100, sd_r = 100, max = 0.9) # sd chose to simulate a width around 450
y.asym_l <- src_asym(x = x, opt = 2000, sd_l = 150, sd_r = 100, max = 0.9) 
y.asym_r <- src_asym(x = x, opt = 2000, sd_l = 100, sd_r = 150, max = 0.9)

# Figure: Representation of the three types of species response curves. 
data.frame("x" = x,
           "y" = c(y.sym, y.asym_l, y.asym_r),
           "toc" = rep(c("sym", "asym_l", "asym_r"), each = length(x))) %>%
  ggplot(aes(x = x, y = y, col = toc)) +
  geom_line(aes(linetype = toc),size = 1.5, alpha = 0.6) +
  ylim(c(0, 1)) +
  scale_color_manual(values = c("red", "darkgreen", "grey20")) +
  scale_linetype_manual(values = c(6, 6, 1)) +
  theme_bw()

## Simulated SRCs for each species in scenario A1 ##
par(mfrow = c(3,3))
for (i in 1:N){
  plot(X_A1[order(X_A1)], psi[,i,1,1], 
       type = "l", 
       xlab = "", ylab = "", 
       ylim = c(0,1.05),
       xlim = range(X_A1),
       col = "blue",
       main = paste(data.sp[i, ]$curve.type1,
                    data.sp[i, ]$curve.type2))
  segments(x0 = data.sp[i, ]$opt,
           x1 = data.sp[i, ]$opt,
           y0 = - 1, 
           y1 = data.sp[i, ]$max.occ1,
           lty = 2,
           col = "blue")
  lines(X_A1[order(X_A1)], psi[,i,2,1],
        col = "red")
  segments(x0 = (data.sp[i, ]$opt + data.sp[i, ]$shift),
           x1 = (data.sp[i, ]$opt + data.sp[i, ]$shift),
           y0 = - 1, 
           y1 = data.sp[i, ]$max.occ2,
           lty = 2,
           col = "red")
  text(x = (data.sp[i, ]$opt + data.sp[i, ]$shift / 2),
       y = 1,
       labels = paste("+", round(data.sp[i, ]$shift, 0), "m"))
  arrows(x0 = (data.sp[i, ]$opt),
         x1 = (data.sp[i, ]$opt + data.sp[i, ]$shift),
         y0 = 0.95, 
         y1 = 0.95,
         length = 0.05)
}

## Simulated SRCs for each species in scenario A1 ##
par(mfrow = c(3,3))
for (i in 1:N){
  plot(X_A2[order(X_A2)], psi[,i,1,2], 
       type = "l", 
       xlab = "", ylab = "", 
       ylim = c(0, 1.05), 
       xlim = range(X_A2),
       col = "blue",
       main = paste(data.sp[i, ]$curve.type1,
                    data.sp[i, ]$curve.type2))
  segments(x0 = data.sp[i, ]$opt,
           x1 = data.sp[i, ]$opt,
           y0 = - 1, 
           y1 = data.sp[i, ]$max.occ1,
           lty = 2,
           col = "blue")
  lines(X_A2[order(X_A2)], psi[,i,2,2],
        col = "red")
  segments(x0 = (data.sp[i, ]$opt + data.sp[i, ]$shift),
           x1 = (data.sp[i, ]$opt + data.sp[i, ]$shift),
           y0 = - 1, 
           y1 = data.sp[i, ]$max.occ2,
           lty = 2,
           col = "red")
  text(x = (data.sp[i, ]$opt + data.sp[i, ]$shift / 2),
       y = 1,
       labels = paste("+", round(data.sp[i, ]$shift, 0), "m"))
  arrows(x0 = (data.sp[i, ]$opt),
         x1 = (data.sp[i, ]$opt + data.sp[i, ]$shift),
         y0 = 0.95, 
         y1 = 0.95,
         length = 0.05)
}

save(list = c("data.sp", "Z", "psi", "X_A1", "X_A2"),
     file = here::here("results", "simulated_data_asymSRC.RData"))
