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
runSim = function(form, data, simPath, simRoot,
                  responseDist = c("gamma", "normal"),
                  genInits = FALSE,
                  nUpdate = 10000, nChains = 1, nIter = 1000){

  if(!dir.exists(simPath)){
    dir.create(file.path(simPath))
  }

  bugsInput = makeBUGSinputs(form, data, responseDist, genInits = genInits)

  ## Write the BUGS file to disk
  bugsFile = file.path(simPath, glue("{simRoot}.bugs.r"))
  writeLines(bugsInput$bugsModelString, bugsFile)


  ## compile the model
  if(genInits){
    system.time({sim = jags.model(file = bugsFile,
                                  data = bugsInput$bugsData,
                                  inits = list(bugsInput$bugsInits),
                                  n.chains = nChains)})
  }else{
    system.time({sim = jags.model(file = bugsFile,
                                  data = bugsInput$bugsData,
                                  n.chains = nChains)})
  }
  update(sim, nUpdate)

  responseDist = match.arg(responseDist)

  parameters = c("pred")

  if(responseDist == "gamma"){
    parameters = c(parameters, "log.Mu", "shape", "rate")
  }else{
    parameters = c(parameters, "Mu")

    if(any(unlist(effects))){ ## any other model than ln-0
      parameters = c(parameters, "mu")
    }
  }

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

  if(effects$bVarEffect){
    parameters = c(parameters, "tau", "tau0")
  }else{
    parameters = c(parameters, "tau")
  }




  sim.sample = coda.samples(sim, parameters, n.iter = nIter)
  simSummary = summary(sim.sample)

  save(sim.sample, simSummary, file = file.path(simPath, glue("{simRoot}.Rda")))
}

