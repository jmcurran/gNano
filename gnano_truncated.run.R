#load required libraries
library(here)
library(rjags)
#loads the make bugs data method
source("makeBugsData.R")
#creates datafile
bugsData = makeBUGSdata()
bugsFile = here("gnano_truncated.bugs.R")

## compile the model
sim = jags.model(file = bugsFile, data = bugsData)

## do a bit of burn in - no idea what is sufficient at this point
system.time(update(sim, 100000))

## The parameters are we interested
parameters = c("pred", "lambda", "mu.dye", "mu.amp", "sigma.sq.dye", "sigma.sq.amp", "T", "loglik")
## run the model
sim.sample = coda.samples(model = sim, variable.names = parameters, n.iter = 50000, thin = 50)

#creates graphs
source("gNano.graphs.R")
saveDir <- here("graphs - uniform//")
createGraphs(sim.sample, saveDir)
