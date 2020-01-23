model
{
  Mu ~ dnorm(0, 1e-05)
  tau ~ dgamma(0.001, 0.001)
  for (i in 1:N) {
    y[i] ~ dlnorm(Mu, tau)
    pred[i] ~ dnorm(Mu, tau)
  }
}
