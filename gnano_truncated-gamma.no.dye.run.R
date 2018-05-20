#must have already run the original model to create the bugsData
bugsFileGamma = here("gnano_truncated-gamma.bugs.R")

## compile the model
simGamma = jags.model(file = bugsFileGamma, data = bugsData)

## do a bit of burn in - no idea what is sufficient at this point
system.time(update(simGamma, 100000))

## What parameters are we interested?
parameters = c("pred", "lambda", "mu.dye", "mu.amp", "sigma.sq.dye", "sigma.sq.amp", "T")
sim.sample.truncated.gamma = coda.samples(model = simGamma, variable.names = parameters, n.iter = 50000, thin = 50)

## now go an print out graphs before doing comparison below
## (need to change calls i ngraphing code from sim.sample to sim.sample.truncated.gamma)

## model comparison code
simGamma = jags.model(file = bugsFileGamma, data = bugsData, n.chain = 4)
sim = jags.model(file = bugsFile, data = bugsData, n.chain = 4)

## now do model comparison using DIC
dic.unif <- dic.samples(model = sim, n.iter = 50000, thin = 50, type = "pD")
dic.gamma <- dic.samples(model = simGamma, n.iter = 50000, thin = 50, type = "pD")
DIC.comparison <- diffdic(dic.unif, dic.gamma)
DIC.comparison


