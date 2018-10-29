buildModel = function(responseDist, bLocusEffect = FALSE, bProfileEffect = FALSE, bDyeEffect = FALSE, bDoseEffect = FALSE){
  locusEffects = "
    #alpha is locus effect
    alpha.mu ~ dnorm(0, 0.000001)
    alpha.tau ~ dgamma(0.001, 0.001)
    alpha.sigma = 1 / sqrt(alpha.tau)
    for(l in 1:(numLoci - 1)){
      alpha.locus[l] ~ dnorm(alpha.mu, alpha.tau)
    }
    alpha.locus[numLoci] = -sum(alpha.locus[1:(numLoci - 1)])
  "

  profileEffects = "
    #beta is profile effect
    beta.mu ~ dnorm(0, 0.000001)
    beta.tau ~ dgamma(0.001, 0.001)
    beta.sigma = 1 / sqrt(beta.tau)
    for(p in 1:(numProfiles - 1)){
      beta.profile[p] ~ dnorm(beta.mu, beta.tau)
    }
    beta.profile[numProfiles] = -sum(beta.profile[1:(numProfiles - 1)])
  "

  dyeEffects = "
    #gamma is profile effect
    gamma.mu ~ dnorm(0, 0.000001)
    gamma.tau ~ dgamma(0.001, 0.001)
    gamma.sigma = 1 / sqrt(gamma.tau)
    for(p in 1:(numDyes - 1)){
      gamma.dye[p] ~ dnorm(gamma.mu, gamma.tau)
    }
    gamma.dye[numProfiles] = -sum(gamma.dye[1:(numDyes - 1)])
  "

  if(responseDist == "gamma"){
    if(!any(c(bLocusEffect, bProfileEffect, bDyeEffect, bDoseEffect))){ #if there are no effects then return the simplest model
      modelString = "model
      {
        log.Mu ~ dnorm(0, 0.000001)
        Mu = exp(log.Mu)
        tau ~ dgamma(0.001, 0.001)
        s = 1/sqrt(tau)
        rate = Mu * tau
        shape = Mu *rate

        for(i in 1:N){
          y[i] ~ dgamma(shape, rate)
        }
      }"
    }else{
      modelString = "model
        {
      "

      meanModel = "y[i] = log.Mu[i]"

      if(bLocusEffect){
        modelString = paste0(modelString, locusEffects)
        meanModel = paste0(meanModel, "+ alpha.locus[locus[i]]")
      }

      if(bProfileEffect){
        modelString = paste0(modelString, profileEffects)
        meanModel = paste0(meanModel, "+ beta.profile[profile[i]]")
      }

      if(bDyeEffect){
        modelString = paste0(modelString, dyeEffects)
        meanModel = paste0(meanModel, "+ gamma[dye[i]]")
      }

      if(bDoseEffect){
        meanModel = paste0(meanModel, "+ log(X[i])")
      }

      main = paste0("
        log.Mu ~ dnorm(0, 0.000001)
        tau ~ dgamma(0.001, 0.001)

        for(i in 1:N){",
        meanModel,
        "
          mu[i] = exp(log.mu[i])
          rate[i] = mu[i] * tau
          shape[i] = mu[i] *rate[i]

          y[i] ~ dgamma(shape[i], rate[i])
          pred[i] ~ dgamma(shape[i], rate[i])
        }")

      modelString = paste0(modelString, main, "\n}")
    }
  }else{
    if(!any(c(bLocusEffect, bProfileEffect, bDyeEffect, bDoseEffect))){ #if there are no effects then return the simplest model
      modelString = "model
      {
          Mu ~ dnorm(0, 0.00001)
          tau ~ dgamma(0.001, 0.001)
          for(i in 1:N){
            log.y[i] ~ dnorm(Mu, tau)
          }"
    }else{
      modelString = "model
        {
      "

      meanModel = "log.y[i] = Mu[i]"

      if(bLocusEffect){
        modelString = paste0(modelString, locusEffects)
        meanModel = paste0(meanModel, "+ alpha.locus[locus[i]]")
      }

      if(bProfileEffect){
        modelString = paste0(modelString, profileEffects)
        meanModel = paste0(meanModel, "+ beta.profile[profile[i]]")
      }

      if(bDyeEffect){
        modelString = paste0(modelString, dyeEffects)
        meanModel = paste0(meanModel, "+ gamma[dye[i]]")
      }

      if(bDoseEffect){
        meanModel = paste0(meanModel, "+ log(X[i])")
      }

      main = paste0("
                    Mu ~ dnorm(0, 0.000001)
                    tau ~ dgamma(0.001, 0.001)

                    for(i in 1:N){",
                    meanModel,
                    "
                    log.y[i] ~ dnorm(mu[i], aph[profile[i]] * tau)
                    pred[i] ~ dnorm(mu[i], aph[profile[i]] * tau)
                    }")

      modelString = paste0(modelString, main, "\n}")
    }
  }

  modelString = paste(tidy_source(text = modelString, arrow = FALSE, indent = 2)$text.tidy, collapse = "\n")
}

#' Make BUGS inputs
#'
#' A simple mechanism for selecting the set of models we want to explore
#'
#' @param form a formula
#' @param data.df a data.frame
#' @param responseDist choice of the peak heights being gamma or (log)-normal
#'
#' @return a list with some stuff in it
#' @export
#'
#' @examples
#' data.df = readData()
#' makeBUGSinputs(data = data.df)
makeBUGSinputs = function(form = formula("y ~ 1"), data.df, responseDist = c("gamma", "normal")){
  if(!is_formula(form)){
    stop("form must be a formula")
  }

  responseDist = match.arg(responseDist)

  vars = attr(terms(form), "term.labels")

  if(!all(grepl("locus|dye|profile|X", vars))){
    mismatch = vars[-grep("locus|dye|profile|X", vars)]
    paste0("The variables : ", mismatch, " don't match the list of allowable variables" )
  }

  ## lazy - always assume y and number of observations are in the reponse
  ## will need to hand logs though
  if(responseDist == "gamma"){
    bugsData = list(y = data.df$obs,
                    N = length(data.df$obs))
  }else{
    aveLogPeakHeight = data.df %>%
      group_by(prof) %>%
      summarise(alph = mean(log(obs), na.rm = TRUE)) %>%
      pull(alph)

    bugsData = list(log.y = log(data.df$obs),
                    N = length(data.df$obs),
                    aveLogPeakHeight = aveLogPeakHeight)
  }

  m = match(c("locus", "profile", "dye", "X"), vars)
  names(m) = c("locus", "profile", "dye", "X")

  if(!is.na(m["locus"])){
    bugsData$locus = data.df$loc
    bugsData$numLoci = 31 ## this may not be true
    bLocusEffect = TRUE
  }else{
    bLocusEffect = FALSE
  }

  if(!is.na(m["profile"])){
    bugsData$profile = data.df$profile
    bugsData$numProfiles = 102
    bProfileEffect = TRUE
  }else{
    bProfileEffect = FALSE
  }

  if(!is.na(m["dye"])){
    bugsData$dye = data.df$dye
    bugsData$numDyes = 4
    bDyeEffect = TRUE
  }else{
    bDyeEffect = FALSE
  }

  if(!is.na(m["X"])){
    bugsData$X = data.df$X
    bDoseEffect = TRUE
  }else{
    bDoseEffect = FALSE
  }

  return(
    list(
      bugsData = bugsData,
      bugsModelString = buildModel(responseDist = responseDist,
                                   bLocusEffect = bLocusEffect,
                                   bProfileEffect = bProfileEffect,
                                   bDyeEffect = bDyeEffect,
                                   bDoseEffect = bDoseEffect),
      responseDist = responseDist,
      modelFormula = form,
      effects = list(bLocusEffect = bLocusEffect,
                     bProfileEffect = bProfileEffect,
                     bDyeEffect = bDyeEffect,
                     bDoseEffect = bDoseEffect)
    )
  )
}
