model{
  for(p in 1:numProfiles){
    mu.profile[p] ~ dnorm(0, 0.00001)
  }
  
  Mu.locus ~ dnorm(0, 0.00001)
  for(l in 1:numLoci){
    mu.locus[l] ~ dnorm(Mu.locus, 0.00001)
  }
  
  for(i in 1:N){
    y[i] ~ dnorm(mu[i], 0.000001)
    mu[i] = mu.profile[profile[i]] + mu.locus[locus[i]]
  }
}