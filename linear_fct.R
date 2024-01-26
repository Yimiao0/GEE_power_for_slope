library(MASS)
library(tidyverse)
library(gee)
library(foreach)
library(doParallel)
# detectCores()


# p is Pr(a subject is in treatment group)
# N is number of subjects
# t is a sequence of time points
# sigma.big is covariance matrix for vector Y_i
data.linear <- function(p, N, t, beta, sigma.big){
  ### assign control/treatment group
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
  colnames(Y) <- c(0, t)
  ### organize the longitudinal data set
  df.wide <- data.frame(as.data.frame(Y),
                        id = 1:N,
                        A = A)
  df.long <- df.wide %>% 
    pivot_longer(cols = X1:X6,
                 names_to = "time", 
                 values_to = "response") %>% 
    rename(Y0 = X0)
  df.long$time <- as.numeric(sub("X", "", df.long$time))
  
  # regression model conditional on Y0
  beta.new <- c(beta, rho) # beta0, beta1, beta2, beta3, beta4=rho
  beta.new[1] <- (1-rho)*beta.new[1] # beta0 = (1-rho)*beta0
  beta.new[1] <- beta.new[1] + rho*mean(Y[,1]) # beta0^* = beta0 + beta4*mean(Y_0)
  
  
  return(list(df.long = df.long,
              beta.ycentered = beta.new,
              A = A,
              Y0 = Y[,1]))
}


# SIMULATION
N <- 20
p <- 0.5 # proportion of subjects in treatment group
t <- c(1, 2, 3, 6) # fixed timepoints
beta <- c(0.5, 1, 1, -0.3) # beta1 + beta3*A as the total slope. beta1:placebo, beta3:treatment prevents disease progression
J <- length(t)
# construct Sigma: covariance matrix
## correlation matrix
rho <- 0.3 # exchangeable correlation for all pairs of measurements within Y_i
corr <- matrix(rho, nrow = J+1, ncol = J+1)
diag(corr) <- 1
## covariance matrix
sigma.small <- 1
sigma.big <- sigma.small^2*corr







# GEE: calculate sample size
df.result <- data.linear(p = p, N = N, t = t, beta = beta, sigma.big = sigma.big)
df.long <- df.result$df.long
beta.ycentered <- df.result$beta.ycentered
A <- df.result$A
Y0 <- df.result$Y0
Y0.centered <- Y0 - mean(Y0)

sigma.small.result <- sigma.small
rho.result <- rho
beta.result <- beta.ycentered # beta0, beta1, beta2, beta3, beta4=rho
D.lst <- list()
for (i in 1:N) {
  D.lst[[i]] <- matrix(NA, nrow = J, ncol = 5, dimnames = list(paste0("t", t), c("intcp", "t", "A", "At", "Y0cent")))
  D.lst[[i]][,1] <- 1
  D.lst[[i]][,2] <- t
  D.lst[[i]][,3] <- A[i]
  D.lst[[i]][,4] <- A[i] * t
  D.lst[[i]][,5] <- Y0.centered[i]
}
R <- matrix(rho.result/(1+rho.result), nrow = J, ncol = J)
diag(R) <- 1
H.lst <- lapply(D.lst, function(x) t(x)%*%solve(R)%*%x)
S <- Reduce("+", H.lst)/N
var.beta <- sigma.small.result^2 * (1 - rho.result^2) / N * solve(S)
var.beta3 <- diag(var.beta)[4]
# var.beta3
# sqrt(var.beta3)
# gee.result <- gee(response ~ time + A + A*time + Y0, id = id, data = df.long,
#                   corstr="exchangeable")
# summary(gee.result)
z.alpha <- qnorm(p = 0.975)
z.b <- qnorm(p = 0.8)
samp.size <- ((1 - rho.result^2)*sigma.small.result^2*(z.alpha + z.b)^2 / beta.result[4]^2) * S[4,4]
samp.size <- ceiling(samp.size)


# calculate power by simulation
numCores <- detectCores()

fx <- function(p, samp.size, t, beta, sigma.big) {
  sim.data.result <- data.linear(p = p, N = samp.size, t = t, 
                                 beta = beta, sigma.big = sigma.big)
  sim.df.long <- sim.data.result$df.long
  
  sim.gee.result <- gee(response ~ time + A + A*time + Y0,
                        id = id, data = sim.df.long,
                        corstr="exchangeable")

  sim.gee.summary <- summary(sim.gee.result)
  sim.est.beta3 <- sim.gee.summary$coefficients["time:A", "Estimate"]
  # sim.var.beta3 <- sim.gee.summary$coefficients["time:A","Robust S.E."]
  
  # hypothesis test
  sim.var.beta3 <- 
  rej <- abs(sim.est.beta3/sqrt(sim.var.beta3)) >= z.alpha

  return(rej)
}


sim.num <- 1000

registerDoParallel(numCores) # use multicore, set to the number of our cores
ptm <- proc.time()
rej <- foreach(icount(sim.num), .combine = cbind) %dopar% {
  fx(p, samp.size, t, beta, sigma.big)
}
tm1 <- proc.time() - ptm

registerDoParallel(numCores) # use multicore, set to the number of our cores
ptm <- proc.time()
rej.compare <- foreach(icount(sim.num), .combine = cbind) %do% {
  fx(p, samp.size, t, beta, sigma.big)
}
tm2 <- proc.time() - ptm

tm1
# user    system  elapsed 
# 49.094  15.398  11.746 
tm2
# user    system  elapsed 
# 28.946  1.585   30.095 

# It is specifically mentioned and illustrated with examples that 
# indeed sometimes it's slower to set this up, because of 
# having to combine the results from the separate parallel processes 
# in the package doParallel.


mean(rej)





##################################################
#### solve beta, sigma.small, rho ################
#### pretend we have solved and obtained them ####
#### Then we're going to calculate Var(beta) #####
#### we'll need (sigma.small, rho, D_i, and R) ###
##################################################
## estimate: sigma.small, beta, rho ##############
## calculate D_i's ###############################

## calculate R ###################################

## calculate H_i = D_iR-1D_i
## calculate S = 1/N*sum(H_i)
## Var(beta)



## calculate sample size N (two-sided type I)
N <- 1000
generate.data <- data.linear(N)
df.long <- generate.data$df.long
beta.ycentered <- generate.data$beta.ycentered

gee.result <- gee(response ~ time + A + A*time + Y0, id = id, data = df.long,
                  corstr="exchangeable")
summary(gee.result)




# ## compare with r package `longpower`:
# library(longpower)
# library(lme4)
# 
# ## definition of pct.change?
# lme.result <- lmer(response ~ time + A + A*time + Y0 + (1 | id), data = df.long) # (1 | id) specifies a random intercept for each subject
# longpower.result <- lmmpower(lme.result, pct.change = 0.3, t = t, power = 0.8, sig.level = 0.05, alternative = "one.sided")
# samp.size.longpower <- longpower.result$N
# samp.size.longpower
# ## ? or
# lmmpower(pct.change = 0.3, t = t, r = rho.result,
#          sigma.small = sigma.small.result,
#          # sig2.i, # estimate of var(random intercept) ?
#          # sig2.s, # estimate of var(random slope)
#          # sig2.e, # estimate of residual variance
#          # cov.s.i, # estimate of covariance of random slope and intercept
#          power = 0.8, sig.level = 0.05, alternative = "one.sided")



## use our samp.size & samp.size.longpower to generate data, then conduct hypothesis test
### *: make the process of generating data a separate function so it can be easily used
### *: make the process of estimating parameters a separate function
## use sample size -> generate data set -> estimate parameters using GEE -> conduct hypothesis test
## -> obtain the test result: reject or not reject (1 or 0)
## repeat the process (from line 118 to 119) 500 times -> calculate power and type I using the 500 simulations

