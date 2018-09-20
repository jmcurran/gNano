#' Run a simulation/fit the model
#'
#' @param form
#' @param data
#' @param simPath
#' @param responseDist
#' @param nUpdate
#' @param nChains
#' @param nIter

#' @export
runSim = function(form, data, simPath, responseDist = c("gamma", "normal"),
                  nUpdate = 10000, nChains = 1, nIter = 1000){

  if(!dir.exists(simPath)){
    dir.create(file.path(simPath))
  }

  bugsInput = makeBUGSinputs(form, data, responseDist)

  ## Write the BUGS file to disk
  bugsFile = file.path(simPath, "model.bugs.r")
  writeLines(bugsInput$bugsModelString, file = bugsFile)


  bugsFile = here("gamma.simple.R")

  ## compile the model
  system.time({sim = jags.model(file = bugsFile,
                                data = bugsInput$bugsData,
                                n.chains = nChains)})
  update(sim, nUpdate)

  parameters = c("pred", "log.Mu", "tau")

  effects = bugsInput$effects
  if(effects$bLocusEffect){
    parameters = c(parameters, "alpha.mu", "alpha.sigma", "alpha.locus")
  }

  if(effects$bProfileEffect){
    parameters = c(parameters, "beta.mu", "beta.sigma", "beta.profile")
  }

  if(effects$bDyeEffect){
    parameters = c(parameters, "gamma.mu", "gamma.sigma", "gamma.dye")
  }


  sim.sample = coda.samples(sim, parameters, n.iter = nIter)
  simSummary = summary(sim.sample)

  save(sim.sample, simSummary, file.path(simPath, "results.rda"))
}
