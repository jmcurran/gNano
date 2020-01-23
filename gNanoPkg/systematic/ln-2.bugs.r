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
  
  Mu ~ dnorm(0, 1e-06)
  tau ~ dgamma(0.001, 0.001)
  
  for (i in 1:N) {
    mu[i] = Mu + beta.profile[profile[i]]
    y[i] ~ dlnorm(mu[i], tau)
    pred[i] ~ dnorm(mu[i], tau)
  }
}
