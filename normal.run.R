#load required libraries
library(here)
library(rjags)

#loads the make bugs data method
source("makeBugsData.R")

#creates datafile
bugsData = makeBUGSdata("K1")
nChains = 1

bugsFile = here("gnano_normal.bugs.R")

## compile the model
system.time({sim = jags.model(file = bugsFile,
                 data = bugsData,
                 n.chains = nChains)})

## do a bit of burn in - no idea what is sufficient at this point
system.time(update(sim, 100000))
#system.time(update(sim, 100))

## The parameters are we interested
parameters = c("pred",
               "lambda",
               "mu.dye",
               "mu.amp",
               "Mu.dye",
               "Mu.amp",
               "T",
               "loglik")

## run the model
system.time({
sim.sample = coda.samples(
  model = sim,
  variable.names = parameters,
  n.iter = 50000,
  thin = 50
)})

#creates graphs
#source("gNano.graphs.R")
#saveDir <- here("graphs - uniform//")
#saveDir <- here(paste("graphs - ", model_descriptor, "//", sep = ""))
#createGraphs(sim.sample, saveDir, bugsData)

#run the WAIC
#source("gNano.WAIC.R")
#calculateWAIC(sim.sample, saveDir, bugsData)

## save sample
#save(sim.sample, file = paste0(saveDir, "sim.sample.Rda"))

#save summary stats
simSummary <- summary(sim.sample)
#pathToFile <- paste(saveDir, "simStatsSummary.txt", sep = "")
#write.table(simSummary[[1]],
#             file = pathToFile,
#             col.names = TRUE,
#             quote = F)
# pathToFile <- paste(saveDir, "simStatsSummaryQuantiles.txt", sep = "")
# write.table(simSummary[[2]],
#             file = pathToFile,
#             col.names = TRUE,
#             quote = F)
