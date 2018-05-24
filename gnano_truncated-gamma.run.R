library(here)
library(rjags)

bugsData = makeBUGSdata()
bugsFileGamma = here("gnano_truncated-gamma.bugs.R")

## compile the model
simGamma = jags.model(file = bugsFileGamma, data = bugsData)

## do a bit of burn in - no idea what is sufficient at this point
system.time(update(simGamma, 100000))

## What parameters are we interested?
parameters = c("pred", "lambda", "mu.dye", "mu.amp", "sigma.sq.dye", "sigma.sq.amp", "T", "loglik")
sim.sample.truncated.gamma = coda.samples(model = simGamma, variable.names = parameters, n.iter = 50000, thin = 50)

## now go an print out graphs before doing comparison below
#gNano.graphs(sim.sample.truncated.gamma, here("GammaGraphs//"))
