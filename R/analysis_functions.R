##############
##
## Functions used in analysis
##
## bastien.mourguiart@gmail.com
## file name: analysis_functions.R
##
###############


#' Optimum estimates with mean comparison method
#'
#' @param z contingence table / presence-absence data of several species
#' @param env vector of environmental values of sampling sites
#'
#' @return mean and confidence interval estimates of species optima
#' @export
#'
#' @examples
mean.method <- function(z, env){
# Compute mean elevation estimates and confidence intervals for each species for one survey
  z[z==0] <- NA
  sp.env <- z*env
  test <- t.test(sp.env)
  res <- c("lwr"=test$conf.int[1], "mean"=test$estimate[[1]], "upr"=test$conf.int[2])
  return(res)
}

#' Mean comparison method
#'
#' @param z contingence table / presence-absence data of several species
#' @param env vector of environmental values of sampling sites
#'
#' @return mean and confidence interval estimates of species shifts
#' @export
#'
#' @examples
mean.method2 <- function(z, env){
# Compute shift in mean elevations estimates and confidence intervals for each species between two surveys
  z[z==0] <- NA

  sp.env <- z*env

  test <- t.test(sp.env[,2], sp.env[,1])

  res <- c("lwr"=test$conf.int[1], "mean"=test$estimate[[1]]-test$estimate[[2]], "upr"=test$conf.int[2])
}



#' Unscale bayesian coefficients
#'
#' @param param Parameter values to be unscale
#' @param x environment vector unscaled
#'
#' @return Unscaled parameter values
#' @export
#'
#' @examples
unscale <- function(param, x) {
  m.x <- mean(x)
  sd.x <- sd(x)

  param*sd.x + m.x
}


#' Convergence indicator for MCMC outputs
#'
#' @param x MCMC output
#'
#' @return Return boolean indicating if all parameters of the MCMC output had converge
#' @export
#'
#' @examples
converge <- function(x) { # Return boolean indicating if all parameters had converge
  r <- x$Rhat # rhat value for each parameter
  l <- lapply(r, function(y) sum(y > 1.1)) # look if some parameter have rhat > 1.1
  sum(l > 0) < 1 # summarise in a TRUE/FALSE that indicates if all parameters had converged
}



#' If interval contains true value
#'
#' @param lwr lower boundary of an interval
#' @param upr upper boundary of an interval
#' @param true true value
#'
#' @return returns if true is contained between lwr and upr
#' @export
#'
#' @examples
contains <- function(lwr, upr, true){

  lwr < true & upr > true

}


#' Interval Score
#'
#' @param alpha confidence degree of the interval
#' @param lwr lower boundary of the interval
#' @param upr upper boundary of the interval
#' @param true true value
#'
#' @return Calculate interval score from CI boundaries and true value
#' @export
#'
#' @examples
ISfunction <- function(alpha = 0.05, lwr, upr, true){

  if (true < lwr){
    a <- (upr-lwr)+(2/alpha)*(lwr-true)
  }
  if (true > upr){
    a <- (upr-lwr)+(2/alpha)*(true-upr)
  }
  if (true < upr & true > lwr){
    a <- (upr-lwr)
  }
  return(a)

}



#' Variance of an estimate
#'
#' @param x a numeric vector
#'
#' @return
#' @export
#'
#' @examples
var2 <- function(x){
  sqrt(var(x) / (dplyr::n() - 1))
}




#' Average performance metrics
#'
#' @param .data data containing performance metrics
#' @param ... column names to be group by
#'
#' @return average performance metrics by selected groups
#' @export
#'
#' @examples
average_metrics <- function(.data, ...) {
  dplyr::group_by(.data, ...) %>%
    dplyr::summarise("m.RMSE" = mean(RMSE),
              "m.SE" = mean(SE),
              "m.bias" = mean(bias),
              "m.coverage" = mean(coverage),
              "m.IS" = mean(IS),
              "sd.RMSE" = var2(RMSE),
              "sd.SE" = var2(SE),
              "sd.bias" = var2(bias),
              "sd.coverage" = var2(coverage),
              "sd.IS" = var2(IS)
              )
}


#' Species Response Curve
#' Simulate a SRC based on defined ecological parameters
#'
#' @param pm mean maximum probability of occupancy
#' @param pinf lower bound mawimum probability of occupancy
#' @param psup upper bound mawimum probability of occupancy
#' @param wm mean width
#' @param winf lower bound width
#' @param wsup upper bound width
#' @param opt optimum
#' @param x gradient
#'
#' @return vector or matrix of occupancy probabilities that reflects a unimodal species response curve with a maximum probability of pm, a width of wm, an optimum of 2000 on a x gradient
#' @export
#'
#' @examples SRC() ; plot(x = 1000:3000, y = SRC())
SRC <- function(pm = 0.925, pinf = NULL, psup = NULL,
                wm = 750, winf = NULL, wsup = NULL,
                opt = 2000,
                x = 1000:3000) {

  p <- c(pm, pinf, psup)
  w <- c(wm, winf, wsup)

  a <- log(p/(1-p))
  tol <- 0.5*w/sqrt(2*(a - log(0.05/0.95)))
  psi <- matrix(0, nrow = length(x), ncol = length(p))

  for (i in 1:length(p)) {
    psi[,i] <- plogis(a[i] - ((x - opt)^2/(2*tol[i]^2)))
  }

  return(psi)

}


#' Convert tolerance in species response curve width
#'
#' @param tol ecological tolerance
#' @param a ehmos parameter alpha
#' @param threshold occupancy probability at which SRC width is calculated
#'
#' @return Species response curve width
#' @export
#'
#' @examples
width_from_tol <- function(tol = tol, a = alpha, threshold = .05){
  2*tol*sqrt(2*(a - log(threshold/(1-threshold))))
}


#' Return mcmc samples of ecological parameters
#' Ecological parameters width, pmax, optimum and shift are derived
#' from the model parameters tol, alpha, opt and shift
#' @param out output of an EHMOS
#' @param scale.x boolean saying if x was standardized in the model
#' @param x environmental vector used in the model
#'
#' @return matrix containing ecological parameters
#' @export
#'
#' @examples
mcmc_ehmos <- function(out = out, scale.x = TRUE, x = alti) {

  if (scale.x == TRUE) {
    if (is.null(x)) {
      return(print("Need to specify a vector of environmental values x"))
    }

    # Unscaled parameters
    alpha <- out$sims.list$alpha
    tol <- (out$sims.list$tol)*sd(x)
    opt <- unscale(out$sims.list$opt, x)
    shift <- out$sims.list$shift*sd(x)
    # Transform tol and alpha in width and pmax
    width <- width_from_tol(tol, alpha, .05)
    pmax <- plogis(alpha)

  }else {

    alpha <- out$sims.list$alpha
    tol <- out$sims.list$tol
    opt <- out$sims.list$opt
    shift <- out$sims.list$shift
    # Transform tol and alpha in width and pmax
    width <- width_from_tol(tol, alpha, .05)
    pmax <- plogis(alpha)

  }

  l <- list("pmax" = pmax,
            "width" = width,
            "opt" = opt,
            "shift" = shift
            )
  return(l)

}



#' Return mcmc samples of ecological parameters
#' Ecological parameters width, pmax, optimum and shift are derived
#' from the model parameters beta
#' @param out output of an glmm
#' @param scale.x boolean saying if x was standardized in the model
#' @param x unscaled environmental vector
#'
#' @return matrix containing ecological parameters
#' @export
#'
#' @examples
mcmc_glmm <- function(out = out, scale.x = TRUE, x = alti) {

  # Model parameters
  b0 <- out$sims.list$b0
  b1 <- out$sims.list$b1
  b2 <- out$sims.list$b2
  b3 <- out$sims.list$b3
  b4 <- out$sims.list$b4
  b5 <- out$sims.list$b5

  if (scale.x == TRUE) {
    if (is.null(x)) {
      return(print("Need to specify a vector of environmental values x"))
    }

    # Scaled parameters
    alpha1 <- b0 - (b1^2) / (4*b2)
    alpha2 <- (b0 + b3) - ((b1 + b4)^2) / (4*(b2 + b5))
    alpha <- array(data = c(alpha1, alpha2),
                   dim = c(nrow(alpha1),
                           ncol(alpha1),
                           2
                   )
    )
    tol1 <- sqrt(-1 / (2 * b2))
    tol2 <- sqrt(-1 / (2 * (b2 + b5)))
    tol <- array(data = c(tol1, tol2),
                   dim = c(nrow(tol1),
                           ncol(tol1),
                           2
                   )
                 )
    opt1 <- -b1 / (2 * b2)
    opt2 <- -(b1 + b4) / (2 * (b2 + b5))
    # Transform tol and alpha in width and pmax, and unscale parameters
    width <- width_from_tol(tol*sd(x), alpha, .05)
    pmax <- plogis(alpha)
    opt <- unscale(opt1, x = x)
    shift <- (opt2 - opt1)*sd(x)

  }else {

    # Scaled parameters
    alpha1 <- b0 - (b1^2) / (4*b2)
    alpha2 <- (b0 + b3) - ((b1 + b4)^2) / (4*(b2 + b5))
    alpha <- array(data = c(alpha1, alpha2),
                   dim = c(nrow(alpha1),
                           ncol(alpha1),
                           2
                   )
    )
    tol1 <- sqrt(-1 / (2 * b2))
    tol2 <- sqrt(-1 / (2 * (b2 + b5)))
    tol <- array(data = c(tol1, tol2),
                 dim = c(nrow(tol1),
                         ncol(tol1),
                         2
                 )
    )
    opt1 <- -b1 / (2 * b2)
    opt2 <- -(b1 + b3) / (2 * (b2 + b5))
    # Transform tol and alpha in width and pmax, and unscale parameters
    width <- width_from_tol(tol, alpha, .05)
    pmax <- plogis(alpha)
    opt <- opt1
    shift <- opt2 - opt1

  }

  l <- list("pmax" = pmax,
            "width" = width,
            "opt" = opt,
            "shift" = shift
  )
  return(l)

}


#' EHMOS Species Response Curve
#' Transform EHMOS parameters in occupancy probabilities to obtain a SRC
#'
#' @param alpha
#' @param x
#' @param opt
#' @param shift
#' @param tol
#'
#' @return
#' @export
#'
#' @examples
ehmos.src <- function(alpha, x, opt, shift, tol){

  plogis(alpha - 0.5*((x - (opt + shift)) / tol)^2)

}


#' GLMM Species Response Curve
#' Transform GLMM parameters in occupancy probabilities to obtain a SRC
#'
#' @param b0
#' @param b1
#' @param b2
#' @param b3
#' @param b4
#' @param b5
#' @param x
#'
#' @return
#' @export
#'
#' @examples
glmm.src <- function(b0, b1, b2, b3, b4, b5, x){

  plogis(b0 + b3 + (b1 + b4)*x + (b2 + b5)*(x^2))

}


#' Predict SRC from EHMOS output
#' Predict Species Response Curve, i.e. evolution of occupancy probabilities along
#' a specified gradient, based on mcmc samples of parameter estimates from
#' an EHMOS model wrote in JAGS
#'
#' @param out output of EHMOS from jagsUI
#' @param scaled specify if environmental variable was scaled in the model
#' @param x.pred environmental vector on which to make predictions of SRC
#' @param full_matrix boolean indicating if the user wants full probability matrix or just summary
#'
#' @return one dataframe of species-specific occupancy probabilities
#' @export
#'
#' @examples
pred.ehmos <- function(out = out,
                       scaled = TRUE,
                       x.pred = alti,
                       full_matrix = FALSE){

  if (scaled == TRUE) {
    x <- scale(x.pred)
  }else{
    x <- x.pred
  }


  alpha <- out$sims.list$alpha
  tol <- out$sims.list$tol
  opt <- out$sims.list$opt
  shift <- out$sims.list$shift

  # Defined the indices
  s <- dim(alpha)[1] # length of mcmc samples
  n <- dim(alpha)[2] # number of species
  j <- length(x) # number of sites

  # Replicate env vector s times
  mat.x <- matrix(rep(x, s),
                  nrow = s,
                  byrow = TRUE)

  full_psi <- array(data = 0,
                    dim = c(s, j, n, 2))

  for (i in 1:n) {

    full_psi[, , i, 1] <- ehmos.src(alpha = alpha[, i, 1],
                                    x = mat.x,
                                    opt = opt[, i],
                                    shift = 0,
                                    tol = tol[, i, 1])

    full_psi[, , i, 2] <- ehmos.src(alpha = alpha[, i, 2],
                                    x = mat.x,
                                    opt = opt[, i],
                                    shift = shift[, i],
                                    tol = tol[, i, 2])

  }

  psi <- apply(full_psi, c(2,3,4), function(x)quantile(x, c(0.025,0.5,0.975)))


  reshape2::melt(psi, varnames = c("quantile", "site", "species", "sampling_occ")) %>%
    tidyr::spread(key = quantile, value = value) %>%
    dplyr::mutate("alti" = rep(x.pred, each=n*2))  -> occu

  if (full_matrix == TRUE){
    return(list(full_psi, occu))
  }else{
    return(occu)
  }

}


#' Predict SRC from GLMM output
#' Predict Species Response Curve, i.e. evolution of occupancy probabilities along
#' a specified gradient, based on mcmc samples of parameter estimates from
#' an c-GLMM model wrote in JAGS
#'
#' @param out output of c-GLMM from jagsUI
#' @param scaled specify if environmental variable was scaled in the model
#' @param x.pred environmental vector on which to make predictions of SRC
#' @param full_matrix boolean indicating if the user wants full probability matrix or just summary
#'
#' @return one dataframe of species-specific occupancy probabilities
#' @export
#'
#' @examples
pred.glmm <- function(out = out,
                       scaled = TRUE,
                       x.pred = alti,
                       full_matrix = FALSE){

  if (scaled == TRUE) {
    x <- scale(x.pred)
  }else{
    x <- x.pred
  }


  b0 <- out$sims.list$b0
  b1 <- out$sims.list$b1
  b2 <- out$sims.list$b2
  b3 <- out$sims.list$b3
  b4 <- out$sims.list$b4
  b5 <- out$sims.list$b5

  # Defined the indices
  s <- dim(b0)[1] # length of mcmc samples
  n <- dim(b0)[2] # number of species
  j <- length(x) # number of sites

  # Replicate env vector s times
  mat.x <- matrix(rep(x, s),
                  nrow = s,
                  byrow = TRUE)

  full_psi <- array(data = 0,
                    dim = c(s, j, n, 2))

  for (i in 1:n) {

    full_psi[, , i, 1] <- glmm.src(b0[, i],
                                   b1[, i],
                                   b2[, i],
                                   0,
                                   0,
                                   0,
                                   mat.x)

    full_psi[, , i, 2] <- glmm.src(b0[, i],
                                   b1[, i],
                                   b2[, i],
                                   b3[, i],
                                   b4[, i],
                                   b5[, i],
                                   mat.x)

  }

  psi <- apply(full_psi, c(2,3,4), function(x)quantile(x, c(0.025,0.5,0.975)))


  reshape2::melt(psi, varnames = c("quantile", "site", "species", "sampling_occ")) %>%
    tidyr::spread(key = quantile, value = value) %>%
    dplyr::mutate("alti" = rep(x.pred, each=n*2))  -> occu

  if (full_matrix == TRUE){
    return(list(full_psi, occu))
  }else{
    return(occu)
  }

}


#' Summarise mcmc samples of ecological parameters
#'
#' @param X mcmc samples of ecological parameters obtained with mcmc_ehmos() or mcmc_glmm()
#'
#' @return
#' @export
#'
#' @examples
summary_mcmc <- function(X){

  lapply(X,
         FUN = function(x){
           if (length(dim(x)) == 3){
             apply(x,
                   c(2, 3),
                   function(y)
                     quantile(y,
                     c(0.025, 0.5, 0.975),
                     na.rm = TRUE)
                   )
             }else{
               apply(x,
                     2,
                     function(y)
                       quantile(y,
                       c(0.025, 0.5, 0.975),
                       na.rm = TRUE)
                     )
             }
           }
         )

}



#' Return width estimates using t-test method
#'
#' @param Z
#' @param alti
#'
#' @return
#' @export
#'
#' @examples
t.test_width <- function(Z, alti){

  Zb <- Z
  Zb[Zb==0] <- NA
  apply(Zb*alti, c(2,3), function(x)(max(x, na.rm=T)-min(x, na.rm=T)))

}
