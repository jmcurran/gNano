
model {
  for (i in 1:N ) {
    dsn[i] <- ( (2/scale)
                * dnorm( (y[i]-locat)/scale , 0 , 1 )
                * pnorm( skew*(y[i]-locat)/scale , 0 , 1 ) )
    spy[i] <- dsn[i] / C
    ones[i] ~ dbern( spy[i] )
    
    u1[i] ~ dnorm(0, 1)
    u2[i] ~ dnorm(0, 1)
    
    z[i] = ifelse(u2[i] < skew * u1[i], u1[i], -u1[i])
    pred[i] = scale * z[i] + location
  }
  scale ~ dgamma(1.105,0.105)
  locat ~ dnorm(0, 0.001)
  skew ~ dnorm(0, 0.001)
}

