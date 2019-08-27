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
  
  Mu ~ dnorm(0, 1e-06)
  scale ~ dgamma(1.105, 0.105)
  skew ~ dnorm(0, 0.001)
  
  for (i in 1:N) {
    location[i] = Mu + alpha.locus[locus[i]] + beta.profile[profile[i]]
    dsn[i] <- ((2/scale) * dnorm((log.y[i] - location[i])/scale, 0, 1) * pnorm(skew * (log.y[i] - location[i])/scale, 0, 1))
    spy[i] <- dsn[i]/C
    ones[i] ~ dbern(spy[i])
    
    u1[i] ~ dnorm(0, 1)
    u2[i] ~ dnorm(0, 1)
    
    z[i] = ifelse(u2[i] < skew * u1[i], u1[i], -u1[i])
    pred[i] = scale * z[i] + location[i]
    
  }
}
