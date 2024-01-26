library(longpower)
source(file = "./linear_fct.R")



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


# GEE:sample size
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
samp.size


# POWER
numCores <- detectCores()

sim.num <- 10000
## parallel computation
numCores <- detectCores()
registerDoParallel(numCores)
rej.parallel <- foreach(sim = 1:sim.num, .combine = cbind) %dopar% {
  fx(p, samp.size, t, beta, sigma.big, z.alpha, seed = sim)
}
mean(rej.parallel, na.rm = T)
sum(is.na(rej.parallel))



# COMPARISON: longpower
## sample size
### Diggle (Book: Analysis of Longitudinal Data, Page29)
diggle.result <- diggle.linear.power(d = 0.1, t = t, 
                                     sigma2 = sigma.small^2, R = rho, 
                                     alternative = "two.sided", power = 0.80)
diggle.samp.size <- ceiling(diggle.result$N)
diggle.samp.size


### Liu, Liang's method


### lmmpower: create a large sample
# lmmpower(beta = 1, pct.change = -0.1, t = t,
#          sig2.i = 0, sig2.s = 0, sig2.e = 0,
#          cov.s.i = 0, 
#          power = 0.8, alternative = "two.sided")
  
library(lme4)

## definition of pct.change?
lme.result <- lmer(response ~ time + A + A*time + Y0 + (1 | id), data = df.long.Y0center) # (1 | id) specifies a random intercept for each subject
longpower.result <- lmmpower(lme.result, pct.change = 0.3, t = t, power = 0.8, sig.level = 0.05, alternative = "one.sided")
samp.size.longpower <- longpower.result$N
samp.size.longpower

liu.liang.linear.power(delta = 0.5, 
                       u = u, 
                       v = v,
                       sigma2 = sigma2, R = rho, 
                       alternative = "two.sided", power = 0.80)





