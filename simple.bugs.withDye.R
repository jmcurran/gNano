model{
  #alpha is locus effect
  for(l in 1:(numLoci - 1)){
    alpha.locus[l] ~ dnorm(0, 0.00001)
  }
  alpha.locus[numLoci] = -sum(alpha.locus[1:(numLoci - 1)])
  
  #beta is profile effect
  for(p in 1:(numProfiles - 1)){
    beta.profile[p] ~ dnorm(0, 0.00001)
  }
  beta.profile[numProfiles] = -sum(beta.profile[1:(numProfiles - 1)])
  
  #gamma is dye effect
  for(f in 1:(numDyes - 1)){
    gamma.dye[f] ~ dnorm(0, 0.00001)
  }
  gamma.dye[numDyes] = -sum(gamma.dye[1:(numDyes - 1)])
  
  Mu ~ dnorm(0, 0.00001)
  for(i in 1:N){
    y[i] ~ dnorm(mu[i], tau[i])
    tau[i] ~ dgamma(0.001, 0.001)
    mu[i] = Mu + alpha.locus[locus[i]] + beta.profile[profile[i]] + gamma.dye[dye[i]]
  }
}