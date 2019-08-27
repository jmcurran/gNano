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
  
  log.Mu ~ dnorm(0, 0.00001)
  
  for (i in 1:N) {
    location[i] = log.Mu + alpha.locus[locus[i]] + beta.profile[profile[i]] + gamma.dye[dye[i]] + X[i]
    
    dsn[i] <- ( (2/scale)
                * dnorm( (y[i] - location[i]) / scale , 0 , 1 )
                * pnorm( skew * (y[i] - location[i]) / scale , 0 , 1 ) )
    spy[i] <- dsn[i] / C
    ones[i] ~ dbern( spy[i] )
  }
  
  scale ~ dgamma(1.105,0.105)
  skew ~ dnorm(0, 0.001)
}
