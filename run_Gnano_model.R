run_Gnano_model <- function(bugs_Model_R_Filename, model_descriptor) {
  #load required libraries
  library(here)
  library(rjags)
  
  #loads the make bugs data method
  source("makeBugsData.R")
  
  #creates datafile
  bugsData = makeBUGSdata()
  
  #bugsFile = here("gnano_truncated.bugs.R")
  bugsFile = here(bugs_Model_R_Filename)
  
  ## compile the model
  sim = jags.model(file = bugsFile,
                   data = bugsData,
                   n.chains = 4)
  
  ## do a bit of burn in - no idea what is sufficient at this point
  system.time(update(sim, 100000))
  #system.time(update(sim, 100))
  
  ## The parameters are we interested
  parameters = c("pred",
                 "lambda",
                 "mu.dye",
                 "mu.amp",
                 "sigma.sq.dye",
                 "sigma.sq.amp",
                 "T",
                 "loglik")
  
  ## run the model
  sim.sample = coda.samples(
    model = sim,
    variable.names = parameters,
    n.iter = 50000,
    thin = 50
  )
  
  #creates graphs
  source("gNano.graphs.R")
  #saveDir <- here("graphs - uniform//")
  saveDir <- here(paste("graphs - ", model_descriptor, "//", sep = ""))
  createGraphs(sim.sample, saveDir, bugsData)
  
  #run the WAIC
  source("gNano.WAIC.R")
  calculateWAIC(sim.sample, saveDir, bugsData)
  
  #save summary stats
  simSummary <- summary(sim.sample)
  pathToFile <- paste(saveDir, "simStatsSummary.txt", sep = "")
  write.table(simSummary[[1]],
              file = pathToFile,
              col.names = TRUE,
              quote = F)
  pathToFile <- paste(saveDir, "simStatsSummaryQuantiles.txt", sep = "")
  write.table(simSummary[[2]],
              file = pathToFile,
              col.names = TRUE,
              quote = F)
}