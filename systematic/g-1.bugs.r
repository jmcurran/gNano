model
{
  
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
    log.mu[i] = log.Mu + beta.profile[profile[i]]
    mu[i] = exp(log.mu[i])
    rate[i] = mu[i] * tau
    shape[i] = mu[i] * rate[i]
    
    y[i] ~ dgamma(shape[i], rate[i])
    pred[i] ~ dgamma(shape[i], rate[i])
  }
}
