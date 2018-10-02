model{
  #alpha is locus effect
  alpha.mu ~ dnorm(0, 0.000001)
  alpha.tau ~ dgamma(0.001, 0.001)
  alpha.sigma = 1 / sqrt(alpha.tau) 
  for(l in 1:(numLoci - 1)){
    alpha.locus[l] ~ dnorm(alpha.mu, alpha.tau)
  }
  alpha.locus[numLoci] = -sum(alpha.locus[1:(numLoci - 1)])
  
  #beta is profile effect
  beta.mu ~ dnorm(0, 0.000001)
  beta.tau ~ dgamma(0.001, 0.001)
  beta.sigma = 1 / sqrt(beta.tau) 
  for(p in 1:(numProfiles - 1)){
    beta.profile[p] ~ dnorm(beta.mu, beta.tau)
  }
  beta.profile[numProfiles] = -sum(beta.profile[1:(numProfiles - 1)])
  
  #gamma is dye effect
  gamma.mu ~ dnorm(0, 0.000001)
  gamma.tau ~ dgamma(0.001, 0.001)
  gamma.sigma = 1 / sqrt(gamma.tau) 
  for(f in 1:(numDyes - 1)){
    gamma.dye[f] ~ dnorm(gamma.mu, gamma.tau)
  }
  gamma.dye[numDyes] = -sum(gamma.dye[1:(numDyes - 1)])
  tau ~ dgamma(0.001, 0.001)
  #s = 1/sqrt(tau[i])
  
  Mu ~ dnorm(0, 0.00001)
  for(i in 1:N){
    log.mu[i] = Mu + X[i] + alpha.locus[locus[i]] + beta.profile[profile[i]] + gamma.dye[dye[i]]
    mu[i] = exp(log.mu[i])
    template[i] = exp(beta.profile[profile[i]])
    rate[i] = mu[i] * tau * template[i]
    shape[i] = mu[i] * rate[i]
    y[i] ~ dgamma(shape[i], rate[i])

  }
}