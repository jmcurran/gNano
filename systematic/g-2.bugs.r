model
{
  
  # alpha is locus effect
  alpha.mu ~ dnorm(0, 1e-06)
  alpha.tau ~ dgamma(0.001, 0.001)
  alpha.sigma = 1/sqrt(alpha.tau)
  for (l in 1:(numLoci - 1)) {
    alpha.locus[l] ~ dnorm(alpha.mu, alpha.tau)
  }
  alpha.locus[numLoci] = -sum(alpha.locus[1:(numLoci - 1)])
  
  # beta is profile effect
  beta.mu ~ dnorm(0, 1e-06)
  beta.tau ~ dgamma(0.001, 0.001)
  beta.sigma = 1/sqrt(beta.tau)
  for (p in 1:(numProfiles - 1)) {
    beta.profile[p] ~ dnorm(beta.mu, beta.tau)
  }
  beta.profile[numProfiles] = -sum(beta.profile[1:(numProfiles - 1)])
  
  log.Mu ~ dnorm(0, 1e-06)
  tau ~ dgamma(0.001, 0.001)
  
  for (i in 1:N) {
    log.mu[i] = log.Mu + alpha.locus[locus[i]] + beta.profile[profile[i]]
    mu[i] = exp(log.mu[i])
    rate[i] = mu[i] * tau
    shape[i] = mu[i] * rate[i]
    
    y[i] ~ dgamma(shape[i], rate[i])
    pred[i] ~ dgamma(shape[i], rate[i])
  }
}
