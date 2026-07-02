# load packages
library(tidyverse) # version 2.0.0
library(ggpubr) # version 0.6.0
library(msm) # version 1.8.2
library(expm) # version 1.0-0
library(gridExtra) # version 2.3
library(RTMB) # version 1.8
# R version 4.4.2

# load file containing functions
source("ctHMM_functions.R")

# load data (obs_all)
# Make sure your data contains patients with measurements of PCT and lactate as well as
# information about their age and sex.


## Preprocessing
# create design matrix
covsMat <- obs_all %>% 
  group_by(age, sex) %>% 
  summarize()

# match rows from design matrix to IDs
covsMat$cov_group <- 1:nrow(covsMat)
obs_all <- obs_all %>% 
  left_join(., covsMat, by = c("age", "sex"))
covsMat <- covsMat[,-3]

covsMat$age <- as.vector(scale(covsMat$age))

covsMat <- as.matrix(cbind(1, covsMat))


## Model Estimation
# define starting values
par = list(
  delta0 = log(c(0.2, 0.1)), # initial state distribution
  mu0.pct = log(c(1.6, 11.6, 19.7)), # mean of gamma distribution for PCT
  mu0.lac = log(c(0.9, 1.7, 8.6)), # mean of gamma distribution for lactate
  sigma0.pct = log(c(2.1, 17.4, 28.1)), # std. of gamma distribution for PCT
  sigma0.lac = log(c(0.3, 0.5, 5.0)), # std. of gamma distribution for lactate
  beta = matrix(c(-0.1, -0.1, -0.3, -0.7, 0.3, 0.2,
                  0.1, 0.1, -0.1, -0.3, 0.2, -0.1), nrow = 6, byrow = T), # beta coefficients
  mass0 = log(c(0.1, 0.5, 0.2)), # group weights
  locs = matrix(c(-7, -9, -8, -10, -10, -10,
                  -3, -15, -15, -13, -10, -10,
                  -15, -11, -9, -9, -10, -10,
                  -5, -8, -7, -11, -10, -10), nrow = 6, byrow = T) # coefficients for groups of RE
)

# define input parameters
data = list(
  cov_group = obs_all$cov_group,
  infection_id = obs_all$infection_id,
  minutes_after_diagnosis = obs_all$minutes_after_diagnosis,
  pct = obs_all$pct,   
  lactate = obs_all$lactate, 
  N = 3, # number of states
  K = 4, # number of groups of RE
  covsMat = covsMat
)

# model estimation
obj <- MakeADFun(llk_hmm, par, silent = TRUE) # objective function
opt <- optim(obj$par, obj$fn, obj$gr, method = "L-BFGS-B", control = list(maxit = 5000, trace = 3))

# save estimated parameters
mod <- obj$report()
sds <- sdreport(obj)


## Results
# calculate AIC and BIC
opt$value
2*opt$value + 2*53 # AIC
2*opt$value + log(length(obs_all$pct))*53 # BIC


# get parameter estimates
N <- 3 
K <- 4 
mu.pct <- mod$mu.pct 
mu.lac <- mod$mu.lac 
sigma.pct <- mod$sigma.pct 
sigma.lac <- mod$sigma.lac 
delta <- mod$delta 
massnat <- mod$massnat 
beta <- matrix(opt$par[15:26], nrow = N*(N-1), byrow = T) 
locs <- matrix(opt$par[30:53], nrow = N*(N-1)) 

Q_array <- array(0, dim = c(N,N,K,nrow(covsMat))) # transition intensity matrix

for (k in 1:K) {
  
  int_beta <- cbind(locs[,k], beta)
  
  z <- 1:N
  t <- 0
  for (i in z) {
    for (j in z[-i]) {
      t <- t + 1
      Q_array[j, i, k, ] <- exp(covsMat %*% int_beta[t, ]) # fill off-diagonal elements
    }
  }
  
  for (i in 1:nrow(covsMat)) {
    diag(Q_array[, , k, i]) <- -rowSums(Q_array[,,k,i]) # fill diagonal elements so that row sum = 0
  }
}


## Calculation of Confidence Intervals
estParams <- opt$par

# get covariance matrix
cov_matrix <- sds$cov.fixed

# get standard deviations
var <- abs(diag(cov_matrix))
sd <- sqrt(var)

# calculate CIs
CIs <- cbind(estParams - qnorm(0.975) * sd, estParams, 
             estParams + qnorm(0.975) * sd)

# mean and standard deviation
exp(CIs[3:14, ])

# lower and upper CI for initial state probability
delta_lower <- c(1, exp(CIs[1:2, 1]))
delta_lower <- delta_lower / sum(delta_lower)
delta_upper <- c(1, exp(CIs[1:2, 3]))
delta_upper <- delta_upper / sum(delta_upper)

delta_lower
delta
delta_upper

# betas
CIs[15:26, ]

# lower and upper CI for mass
mass_lower <- c(1, exp(CIs[27:29, 1]))
mass_lower <- mass_lower / sum(mass_lower)
mass_upper <- c(1, exp(CIs[27:29, 3]))
mass_upper <- mass_upper / sum(mass_upper)

mass_lower
massnat
mass_upper

# random effects
CIs[30:53,]


## Local State Decoding
# calculate individual weights for random effects
ind_weights_all <- ind_weight_RE(mod, opt, data)

# state decoding
stateProbs <- state_probs_hmm(mod, data, ind_weights_all)

ild <- rep(NA, nrow(obs_all))
for (i in 1:nrow(obs_all)) {
  ild[i] <- which.max(stateProbs[i,])
}
obs_all$localStates <- ild
obs_all <- cbind(obs_all, stateProbs)
colnames(obs_all) <- c("infection_id", "pct", "lactate", "age", "sex", "minutes_after_diagnosis", "cov_group",
                       "local_states", "local_prob_state1", "local_prob_state2", "local_prob_state3")
