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
  
  # gamma is profile effect
  gamma.mu ~ dnorm(0, 1e-06)
  gamma.tau ~ dgamma(0.001, 0.001)
  gamma.sigma = 1/sqrt(gamma.tau)
  for (p in 1:(numDyes - 1)) {
    gamma.dye[p] ~ dnorm(gamma.mu, gamma.tau)
  }
  gamma.dye[numDyes] = -sum(gamma.dye[1:(numDyes - 1)])
  # delta is the locus-dye interaction
  delta.mu ~ dnorm(0, 1e-06)
  delta.tau ~ dgamma(0.001, 0.001)
  delta.sigma = 1/sqrt(delta.tau)
  for (p in 1:(numLocusDyes - 1)) {
    delta.locus.dye[p] ~ dnorm(delta.mu, delta.tau)
  }
  
  delta.locus.dye[numLocusDyes] = -sum(delta.locus.dye[1:(numLocusDyes - 1)])
  
  log.Mu ~ dnorm(0, 1e-06)
  tau ~ dgamma(0.001, 0.001)
  
  for (i in 1:N) {
    log.mu[i] = log.Mu + alpha.locus[locus[i]] + beta.profile[profile[i]] + gamma.dye[dye[i]] + delta.locus.dye[locus.dye[i]] + X[i]
    mu[i] = exp(log.mu[i])
    rate[i] = mu[i] * tau
    shape[i] = mu[i] * rate[i]
    
    y[i] ~ dgamma(shape[i], rate[i])
    pred[i] ~ dgamma(shape[i], rate[i])
  }
}
