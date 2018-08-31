model{
  log.mu ~ dnorm(0, 0.000001)
  mu = exp(log.mu)
  tau ~ dgamma(0.001, 0.001)
  for(i in 1:N){
    y[i] ~
    
  }
}