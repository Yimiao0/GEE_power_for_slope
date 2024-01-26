library(MASS)
library(tidyverse)
library(gee)
library(doParallel)


# FUNCTIONS
## generate data
data.linear <- function(p, N, t, beta, sigma.big){
  # p is Pr(a subject is in treatment group)
  # N is number of subjects
  # t is a sequence of time points
  # sigma.big is covariance matrix for vector Y_i
  
  # assign control/treatment group
  treat.group <- sample(1:N, round(p*N), replace = F)
  A <- rep(0, N)
  A[treat.group] <- 1
  
  J <- length(t)
  
  mu0 <- rep(beta[1], N) # mu_i0 = beta0
  mu <- matrix(NA, nrow = N, ncol = J) # mu_ij, j = 1,...,J
  for (j in 1:J) {
    X <- matrix(1, nrow = N, ncol = 4) # 1, tj, Ai, Aitj
    X[,2] <- t[j]
    X[,3] <- A
    X[,4] <- A*t[j]
    
    mu[,j] <- X%*%beta
  }
  
  # Y matrix: N*(1+J)
  residual <- mvrnorm(n = N, mu = rep(0, J+1), Sigma = sigma.big)
  Y <- cbind(mu0, mu) + residual
  colnames(Y) <- sapply(c(0, t), function(x) paste("Y", x, sep = ""))
  
  # center Y0
  Y0.centered <- Y[,1] - mean(Y[,1])
  
  # organize the longitudinal data set given Y0
  df.wide <- data.frame(id = 1:N,
                        A = A,
                        as.data.frame(Y))
  df.long <- df.wide %>% 
    pivot_longer(cols = Y1:Y6,
                 names_to = "time", 
                 values_to = "response")
  df.long$time <- as.numeric(sub("Y", "", df.long$time))
  
  # df.long with centered Y0
  df.long.Y0center <- df.long
  df.long.Y0center$Y0 <- rep(Y0.centered, each = J)
  
  # regression model conditional on Y0: so coefficients have changed
  beta.new <- c(beta, rho) # beta4=rho
  beta.new[1] <- (1-rho)*beta.new[1] # beta0 = (1-rho)*beta0
  beta.new[1] <- beta.new[1] + rho*mean(Y[,1]) # beta0^* = beta0 + beta4*mean(Y_0)
  names(beta.new) <- c("intercept", "t", "A", "At", "Y_center")
  
  
  return(list(df.long = df.long,
              df.long.Y0center = df.long.Y0center,
              beta.ycentered = beta.new,
              A = A,
              Y0 = Y[,1]))
}

## calculate S matrix for variance(beta)
calc.S <- function(t, A, rho.result, Y0.centered) {
  J <- length(t)
  N <- length(A)
  R <- matrix(rho.result / (1 + rho.result), nrow = J, ncol = J)
  diag(R) <- 1
  D.lst <- list()
  for (i in 1:N) {
    D.lst[[i]] <- matrix(NA, nrow = J, ncol = 5, 
                         dimnames = list(paste0("t", t), 
                                         c("intcp", "t", "A", "At", "Y0cent")))
    D.lst[[i]][,1] <- 1
    D.lst[[i]][,2] <- t
    D.lst[[i]][,3] <- A[i]
    D.lst[[i]][,4] <- A[i] * t
    D.lst[[i]][,5] <- Y0.centered[i]
  }
  H.lst <- lapply(D.lst, function(x) t(x)%*%solve(R)%*%x)
  S <- Reduce("+", H.lst) / N
  
  return(S)
}



# POWER

fx <- function(p, samp.size, t, beta, sigma.big, z.alpha, seed) {
  set.seed(seed)
  J <- length(t)
  sim.df.result <- data.linear(p = p, N = samp.size, t = t, 
                               beta = beta, sigma.big = sigma.big)
  sim.df.long.Y0center <- sim.df.result$df.long.Y0center
  sim.A <- sim.df.result$A
  sim.Y0 <- sim.df.result$Y0
  sim.Y0.centered <- sim.Y0 - mean(sim.Y0)
  
  sim.gee.result <- gee(response ~ time + A + A*time + Y0,
                        id = id, data = sim.df.long.Y0center,
                        corstr="exchangeable")
  sim.residuals <- sim.gee.result$residuals
  
  ######### NOTICE: here sigma.small and rho are conditional values given Y0
  ######### that is, sigma.small = sigma.small^2 * (1-rho^2), rho = rho/(1+rho)
  sim.gee.summary <- summary(sim.gee.result)
  rej <- abs(sim.gee.summary$coefficients["time:A", "Naive z"]) >= z.alpha
  
  return(rej)
}









