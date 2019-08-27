model
{
  log.Mu ~ dnorm(0, 1e-06)
  Mu = exp(log.Mu)
  tau ~ dgamma(0.001, 0.001)
  s = 1/sqrt(tau)
  rate = Mu * tau
  shape = Mu * rate
  
  for (i in 1:N) {
    y[i] ~ dgamma(shape, rate)
    pred[i] ~ dgamma(shape, rate)
  }
}
