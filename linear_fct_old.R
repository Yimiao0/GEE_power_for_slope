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




# PARAMETER SETTING
p <- 0.5
N <- 20
t <- c(1, 2, 3, 6)
beta <- c(0.5, 1, 1, -0.1) # beta1 + beta3*A as the total slope. beta1:placebo, beta3:treatment prevents disease progression
J <- length(t)
## construct Sigma: covariance matrix
rho <- 0.5 # exchangeable correlation for all pairs of measurements within Y_i
corr <- matrix(rho, nrow = J+1, ncol = J+1)
diag(corr) <- 1
sigma.small <- 2
sigma.big <- sigma.small^2 * corr


# GENERATE DATA SET
df.result <- data.linear(p = p, N = N, t = t, beta = beta, sigma.big = sigma.big)
df.long <- df.result$df.long
df.long.Y0center <- df.result$df.long.Y0center
beta.ycentered <- df.result$beta.ycentered
A <- df.result$A
Y0 <- df.result$Y0
Y0.centered <- Y0 - mean(Y0)


# GEE
## for now, use the true value of parameters to calculate sample size
z.alpha <- qnorm(p = 0.975)
z.b <- qnorm(p = 0.8)

beta.result <- beta.ycentered
delta <- beta.result[4]
sigma.small.result <- sigma.small
rho.result <- rho
S <- calc.S(t, A, rho.result, Y0.centered)

samp.size <- ((1 - rho.result^2) * sigma.small.result^2 * (z.alpha + z.b)^2 / delta^2) * solve(S)[4,4]
samp.size <- unname(ceiling(samp.size))


# POWER
numCores <- detectCores()

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
  
  # ######### NOTICE: here sigma.small and rho are conditional value given Y0
  # ######### that is, sigma.small = sigma.small^2 * (1-rho^2), rho = rho/(1+rho)
  # # sigma.small^2 | Y0
  # sim.conditional.sigma.small2.est <- sum(sim.residuals^2) / (samp.size*J - 4)
  # ## or use result given by gee() 
  # # sim.conditional.sigma.small2.est <- sim.gee.result$scale
  # # rho | Y0
  # sim.residuals.mat <- matrix(sim.residuals, nrow = J)
  # multiplied.res.mat <- matrix(0, nrow = J*(J - 1)/2, ncol = ncol(sim.residuals.mat))
  # row.ind <- 0
  # for (prev in 1:(J - 1)) {
  #   for (latter in (prev+1):J) {
  #     row.ind <- row.ind + 1
  #     multiplied.res.mat[row.ind,] <- sim.residuals.mat[prev,] * sim.residuals.mat[latter,]
  #   }
  # }
  # sim.conditional.rho <- sum(multiplied.res.mat) / sim.conditional.sigma.small2.est / (samp.size*J*(J-1)/2 - 4)
  # ## or use result given by gee() 
  # # sim.working.correlation <- sim.gee.result$working.correlation
  # # sim.conditional.rho <- sim.working.correlation[1,2]
  # 
  # ## Q to solve: beta3=-1, set.seed(12), sim.conditional.rho is greater than 0.5 --> sim.rho.est greater than 1
  # sim.rho.est <- sim.conditional.rho / (1 - sim.conditional.rho)
  # sim.sigma.small2.est <- sim.conditional.sigma.small2.est / (1 - sim.rho.est^2)
  # 
  # sim.gee.summary <- summary(sim.gee.result)
  # sim.est.beta3 <- sim.gee.summary$coefficients["time:A", "Estimate"]
  # # sim.var.beta3 <- sim.gee.summary$coefficients["time:A","Robust S.E."]
  # 
  # # Hypothesis test
  # sim.S <- calc.S(t, sim.A, sim.rho.est, sim.Y0.centered)
  # sim.var.beta <- (sim.sigma.small2.est * (1 - sim.rho.est^2) / samp.size) * solve(sim.S)
  # sim.var.beta3 <- sim.var.beta[4,4]
  # rej <- abs(sim.est.beta3/sqrt(sim.var.beta3)) >= z.alpha
  
  sim.gee.summary <- summary(sim.gee.result)
  rej <- abs(sim.gee.summary$coefficients["time:A", "Naive z"]) >= z.alpha
  
  return(rej)
}


sim.num <- 10000
# for loop: takes long time because gee prints too much and cannot be suppressed
ptm1 <- proc.time()

rej.for_loop <- c()
for (sim in 1:sim.num) {
  rej.for_loop[sim] <- fx(p, samp.size, t, beta, sigma.big, z.alpha, seed = sim)
}
mean(rej.for_loop, na.rm = T)
sum(is.na(rej.for_loop))

tm1 <- proc.time() - ptm1


# parallel computation
ptm2 <- proc.time()

numCores <- detectCores()
registerDoParallel(numCores)
rej.parallel <- foreach(sim = 1:sim.num, .combine = cbind) %dopar% {
  fx(p, samp.size, t, beta, sigma.big, z.alpha, seed = sim)
}
mean(rej.parallel, na.rm = T)
sum(is.na(rej.parallel))

tm2 <- proc.time() - ptm2

# compare the time
tm1
tm2




