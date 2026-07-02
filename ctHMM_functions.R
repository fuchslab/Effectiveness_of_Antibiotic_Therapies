# load packages
library(expm) # version 1.0-0
library(LaMa)# version 2.0.6
# R version 4.4.2


## Log-likelihood Function
llk_hmm <- function(par){
  getAll(par, data)
  mu.pct <- exp(mu0.pct); REPORT(mu.pct); ADREPORT(mu.pct)
  mu.lac <- exp(mu0.lac); REPORT(mu.lac); ADREPORT(mu.lac)
  sigma.pct <- exp(sigma0.pct); REPORT(sigma.pct); ADREPORT(sigma.pct)
  sigma.lac <- exp(sigma0.lac); REPORT(sigma.lac); ADREPORT(sigma.lac)
  delta <- c(1, exp(delta0))
  delta <- delta/sum(delta); REPORT(delta); ADREPORT(delta)
  massnat <- c(1, exp(mass0))
  massnat <- massnat/sum(massnat); REPORT(massnat); ADREPORT(massnat)
  
  Q_array <- array(0, dim = c(N,N,K,nrow(covsMat))) # create 4-dim array that contains all Q matrices for each random effect and for each group of age and sex
  
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
  
  # loop over each patient
  mllk.all <- 0
  patient <- unique(infection_id) # get unique IDs
  for (i in 1:length(patient)) { # for each unique ID
    idx <- which(infection_id == patient[i]) # filter data of this ID
    n <- length(idx) # get number of rows of ID
    timeDiff <- diff(minutes_after_diagnosis[idx])
    
    # P matrix
    one.patient<- matrix(rep(NA,K*2),ncol=2)
    allprobs <- matrix(1, nrow = n, ncol = N) # create matrix with N columns and n rows
    ind.pct <- which(!is.na(pct[idx]))
    ind.lac <- which(!is.na(lactate[idx]))
    for (j in 1:N){
      pct.prob <- lac.prob <- rep(1, n)
      pct.prob[ind.pct] <- dgamma(pct[idx[ind.pct]], shape = mu.pct[j]^2/(sigma.pct[j])^2, scale = sigma.pct[j]^2/mu.pct[j])
      lac.prob[ind.lac] <- dgamma(lactate[idx[ind.lac]], shape = mu.lac[j]^2/(sigma.lac[j])^2, scale = sigma.lac[j]^2/mu.lac[j])
      allprobs[,j] <- pct.prob*lac.prob
    }
    
    for (k in 1:K) { # loop over number of random effects
      
      Q <- Q_array[,, k, unique(cov_group[idx])] # get Q matrix for each group of age and sex
      
      Gamma <- tpm_cont(Q, timeDiff)
      
      # forward algorithm
      lscale <- forward_g(delta, Gamma, allprobs)
      
      one.patient[k,1]<- lscale 
      one.patient[k,2]<- massnat[k]
    }
    ma <- max(one.patient[,1]) # necessary to avoid numerical underflow   
    one.patient[,1] <- exp(one.patient[,1]-ma)  
    mllk.all <- mllk.all+(log(sum(one.patient[,1]*one.patient[,2]))+ma)
  }
  return(-mllk.all)
}




## Individual-specific Weights for RE
ind_weight_RE <- function(mod, opt, data){
  getAll(data)
  mu.pct <- mod$mu.pct 
  mu.lac <- mod$mu.lac 
  sigma.pct <- mod$sigma.pct 
  sigma.lac <- mod$sigma.lac 
  delta <- mod$delta 
  massnat <- mod$massnat 
  beta <- matrix(opt$par[15:26], nrow = N*(N-1), byrow = T) 
  locs <- matrix(opt$par[30:53], nrow = N*(N-1)) 
  
  Q_array <- array(0, dim = c(N,N,K,nrow(covsMat)))
  
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
      diag(Q_array[, , k, i]) <- -rowSums(Q_array[,,k,i]) # fill diagonal elements so that rowsum = 0
    }
  }
  
  # loop over each patient
  patient <- unique(infection_id) # get unique IDs
  ind_weights_all <- matrix(NA, nrow = length(patient), ncol = K)
  for (i in 1:length(patient)) { # for each unique ID
    idx <- which(infection_id == patient[i]) # filter data of this ID
    n <- length(idx) # get number of rows of ID
    timeDiff <- diff(minutes_after_diagnosis[idx])
    
    # P matrix
    one.patient<- matrix(rep(NA,K*2),ncol=2)
    allprobs <- matrix(1, nrow = n, ncol = N) # create matrix with N columns and n rows
    ind.pct <- which(!is.na(pct[idx]))
    ind.lac <- which(!is.na(lactate[idx]))
    for (j in 1:N){
      pct.prob <- lac.prob <- rep(1, n)
      pct.prob[ind.pct] <- dgamma(pct[idx[ind.pct]], shape = mu.pct[j]^2/(sigma.pct[j])^2, scale = sigma.pct[j]^2/mu.pct[j])
      lac.prob[ind.lac] <- dgamma(lactate[idx[ind.lac]], shape = mu.lac[j]^2/(sigma.lac[j])^2, scale = sigma.lac[j]^2/mu.lac[j])
      allprobs[,j] <- pct.prob*lac.prob
    }
    
    for (k in 1:K) { # loop over number of random effects
      
      Q <- Q_array[,, k, unique(cov_group[idx])] # get Q matrix for each group of age and sex
      
      Gamma <- tpm_cont(Q, timeDiff)
      
      # forward algorithm
      lscale <- forward_g(delta, Gamma, allprobs)
      
      one.patient[k,1]<- lscale 
      one.patient[k,2]<- massnat[k]
    }
    
    w_llk <- one.patient[,1] + log(one.patient[,2]) # llk already in log scale, weights in prob scale, therefore log
    ma <- max(w_llk) # necessary to avoid numerical underflow   
    w_L <- exp(w_llk - ma)/sum(exp(w_llk - ma))
    
    ind_weights_all[i,] <- w_L
  }
  return(ind_weights_all)
}




## Local Decoding via Forward-backward Algorithm
# forward probabilities
lforward_hmm <- function(mod, n, N, delta, allprobs, Q, timeDiff){
  n             <- n
  lalpha        <- matrix(NA, N, n)
  foo           <- delta * allprobs[1, ]
  sumfoo        <- sum(foo)
  lscale        <- log(sumfoo)
  foo           <- foo / sumfoo
  lalpha[, 1]   <- lscale + log(foo)
  
  for (i in 2:n) {
    Gamma       <- expm(Q * timeDiff[i-1])
    foo         <- foo %*% Gamma * allprobs[i, ]
    sumfoo      <- sum(foo)
    lscale      <- lscale + log(sumfoo)
    foo         <- foo / sumfoo
    lalpha[, i] <- log(foo) + lscale
  }
  return(lalpha)
}


# backward probabilities
lbackward_hmm <- function(mod, n, N, allprobs, Q, timeDiff){
  n             <- n
  lbeta         <- matrix(NA, N, n)
  lbeta[, n]    <- rep(0, N)
  foo           <- rep(1/N, N)
  lscale        <- log(N)
  
  for (i in (n-1):1) {
    Gamma       <- expm(Q * timeDiff[i])
    foo         <- Gamma %*% (allprobs[i + 1, ] * foo)
    lbeta[, i]  <- log(foo) + lscale
    sumfoo      <- sum(foo)
    foo         <- foo/sumfoo
    lscale      <- lscale + log(sumfoo)
  }
  return(lbeta) 
}


# state probabilities
state_probs_hmm <- function(mod, data, ind_weights_all){
  getAll(data)
  mu.pct <- mod$mu.pct
  mu.lac <- mod$mu.lac
  sigma.pct <- mod$sigma.pct
  sigma.lac <- mod$sigma.lac
  delta <- mod$delta
  massnat <- mod$massnat
  beta <- matrix(opt$par[15:26], nrow = N*(N-1), byrow = T)
  locs <- matrix(opt$par[30:53], nrow = N*(N-1))
  
  stateProbs_all <- matrix(NA, ncol = N)
  
  Q_array <- array(0, dim = c(N,N,K,nrow(covsMat)))
  
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
      diag(Q_array[, , k, i]) <- -rowSums(Q_array[,,k,i]) # fill diagonal elements so that rowsum = 0
    }
  }
  
  
  # loop over each patient
  patient <- unique(infection_id) # get unique IDs
  for (i in 1:length(patient)) { # for each unique ID
    idx <- which(infection_id == patient[i]) # filter data of this ID
    n <- length(idx) # get number of rows of ID
    timeDiff <- diff(minutes_after_diagnosis[idx])
    
    # P matrix
    one.patient<- matrix(rep(NA,K*2),ncol=2)
    allprobs <- matrix(1, nrow = n, ncol = N) # create matrix with N columns and n rows
    ind.pct <- which(!is.na(pct[idx]))
    ind.lac <- which(!is.na(lactate[idx]))
    for (j in 1:N){
      pct.prob <- lac.prob <- rep(1, n)
      pct.prob[ind.pct] <- dgamma(pct[idx[ind.pct]], shape = mu.pct[j]^2/(sigma.pct[j])^2, scale = sigma.pct[j]^2/mu.pct[j])
      lac.prob[ind.lac] <- dgamma(lactate[idx[ind.lac]], shape = mu.lac[j]^2/(sigma.lac[j])^2, scale = sigma.lac[j]^2/mu.lac[j])
      allprobs[,j] <- pct.prob*lac.prob
    }
    
    state_probs <- array(NA, dim = c(n,N,K))
    for (k in 1:K) { # loop over number of random effects
      
      Q <- Q_array[,, k, unique(cov_group[idx])] # get Q matrix for each group of age and sex
      
      la <- lforward_hmm(mod, n, N, delta, allprobs, Q, timeDiff)
      lb <- lbackward_hmm(mod, n, N, allprobs, Q, timeDiff)
      
      c <- max(la[, n])
      llk <- c + log(sum(exp(la[, n] - c)))
      for (t in 1:n) {
        state_probs[t,,k] <- exp(la[, t] + lb[, t] - llk)
      }
    }
    weighted_state_probs <- matrix(NA, ncol = N, nrow = n)
    for (j in 1:n) {
      weighted_state_probs[j,] <- t(state_probs[j,,] %*% ind_weights_all[i,])
    }
    stateProbs_all <- rbind(stateProbs_all, weighted_state_probs)
  }
  return(stateProbs_all[-1, ])
}
