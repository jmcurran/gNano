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
  log.Mu ~ dnorm(0, 0.00001)
  
  for (i in 1:N) {
    location[i] = log.Mu + beta.profile[profile[i]]
    
    dsn[i] <- ( (2/scale)
                * dnorm( (y[i] - location[i]) / scale , 0 , 1 )
                * pnorm( skew * (y[i] - location[i]) / scale , 0 , 1 ) )
    spy[i] <- dsn[i] / C
    ones[i] ~ dbern( spy[i] )
  }
  
  scale ~ dgamma(1.105,0.105)
  skew ~ dnorm(0, 0.001)
}
