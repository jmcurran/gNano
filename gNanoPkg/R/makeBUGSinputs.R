buildModel = function(responseDist, bLocusEffect = FALSE, bProfileEffect = FALSE, bDyeEffect = FALSE, bDoseEffect = FALSE,
                      bVarEffect = FALSE){
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
    gamma.dye[numDyes] = -sum(gamma.dye[1:(numDyes - 1)])
  "

  if(responseDist == "gamma"){
    if(!any(c(bLocusEffect, bProfileEffect, bDyeEffect, bDoseEffect))){ #if there are no effects then return the simplest model
      if(!!bVarEffect){
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
            pred[i] ~ dgamma(shape, rate)
          }
        }"
      }else{
        modelString = "model
        {
          log.Mu ~ dnorm(0, 0.000001)
          Mu = exp(log.Mu)
          tau0 ~ dgamma(0.001, 0.001)
          tau1 ~ dgamma(0.001, 0.001)

          for(i in 1:N){
            y[i] ~ dgamma(shape[i], rate[i])
            pred[i] ~ dgamma(shape[i], rate[i])

            tau[i] = aph[profile[i]] / (aph[profile[i]] / tau0 + 1 / tau1)

            rate[i] = Mu * tau[i]
            shape[i] = Mu * rate[i]
          }
        }"
      }
    }else{
      modelString = "model
        {
      "

      meanModel = "log.mu[i] = log.Mu"

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
        meanModel = paste0(meanModel, "+ gamma.dye[dye[i]]")
      }

      if(bDoseEffect){
        meanModel = paste0(meanModel, "+ X[i]")
      }

      main = if(!bVarEffect){
        paste0("
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
      }else{
        paste0("
          log.Mu ~ dnorm(0, 0.000001)
          tau0 ~ dgamma(0.001, 0.001)
          tau1 ~ dgamma(0.001, 0.001)

          for(i in 1:N){",
               meanModel,
               "
            mu[i] = exp(log.mu[i])

            tau[i] = aph[profile[i]] / (aph[profile[i]] / tau0 + 1 / tau1)
            rate[i] = mu[i] * tau[i]
            shape[i] = mu[i] *rate[i]

            y[i] ~ dgamma(shape[i], rate[i])
            pred[i] ~ dgamma(shape[i], rate[i])
          }")
      }

      modelString = paste0(modelString, main, "\n}")
    }
  }else{
    if(!any(c(bLocusEffect, bProfileEffect, bDyeEffect, bDoseEffect))){ #if there are no effects then return the simplest model
      if(!bVarEffect){
        modelString = "model
        {
            Mu ~ dnorm(0, 0.00001)
            tau ~ dgamma(0.001, 0.001)
            for(i in 1:N){
              log.y[i] ~ dnorm(Mu, tau)
              pred[i] ~ dnorm(Mu, tau)
            }
        }"
      }else{
        modelString = "model
        {
            Mu ~ dnorm(0, 0.00001)
            tau0 ~ dgamma(0.001, 0.001)
            tau1 ~ dgamma(0.001, 0.001)

            for(i in 1:N){
              log.y[i] ~ dnorm(Mu, tau[i])

              tau[i] = aph[profile[i]] / (aph[profile[i]] / tau0 + 1 / tau1)

              pred[i] ~ dnorm(Mu, tau[i])
            }
        }"
      }
    }else{
      modelString = "model
        {
      "

      meanModel = "mu[i] = Mu"

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
        meanModel = paste0(meanModel, "+ gamma.dye[dye[i]]")
      }

      if(bDoseEffect){
        meanModel = paste0(meanModel, "+ X[i]")
      }

      main = if(!bVarEffect){
        paste0("
                      Mu ~ dnorm(0, 0.000001)
                      tau ~ dgamma(0.001, 0.001)

                      for(i in 1:N){",
                      meanModel,
                      "
                      log.y[i] ~ dnorm(mu[i], tau)
                      pred[i] ~ dnorm(mu[i], tau)
                      }")
        }else{
          paste0("
                      Mu ~ dnorm(0, 0.000001)
                      tau0 ~ dgamma(0.001, 0.001)
                      tau1~ dgamma(0.001, 0.001)

                      for(i in 1:N){",
                 meanModel,
                 "
                      log.y[i] ~ dnorm(mu[i], tau[i])

                      tau[i] = aph[profile[i]] / (aph[profile[i]] / tau0 + 1 / tau1)

                      pred[i] ~ dnorm(mu[i], tau[i])
                      }")
        }
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
makeBUGSinputs = function(form = formula("y ~ 1"), data.df, responseDist = c("gamma", "normal"),
                          genInits = FALSE){
  if(!is_formula(form)){
    stop("form must be a formula")
  }

  responseDist = match.arg(responseDist)

  vars = attr(terms(form), "term.labels")

  if(!all(grepl("locus|dye|profile|X|V", vars))){
    mismatch = vars[-grep("locus|dye|profile|X|V", vars)]
    paste0("The variables : ", mismatch, " don't match the list of allowable variables" )
  }

  ## lazy - always assume y and number of observations are in the reponse
  ## will need to handle logs though
  if(responseDist == "gamma"){
    aveLogPeakHeight = data.df %>%
      group_by(prof) %>%
      summarise(alph = mean(log(obs), na.rm = TRUE)) %>%
      pull(alph)

    bugsData = list(y = data.df$obs,
                    N = length(data.df$obs),
                    aph = aveLogPeakHeight)
  }else{
    aveLogPeakHeight = data.df %>%
      group_by(prof) %>%
      summarise(alph = mean(log(obs), na.rm = TRUE)) %>%
      pull(alph)

    bugsData = list(log.y = log(data.df$obs),
                    N = length(data.df$obs),
                    aph = aveLogPeakHeight)
  }

  m = match(c("locus", "profile", "dye", "X", "V"), vars)
  names(m) = c("locus", "profile", "dye", "X", "V")

  bugsInits = list()

  if(!is.na(m["locus"])){
    bugsData$locus = data.df$loc
    bugsData$numLoci = 31 ## this may not be true
    bLocusEffect = TRUE

    if(genInits){
      bugsInits$alpha.mu = 0
      bugsInits$alpha.tau = 1
      #bugsInits$alpha.locus = rep(0, bugsData$numLoci - 1)
    }

  }else{
    bLocusEffect = FALSE
  }

  if(!is.na(m["profile"])){
    bugsData$profile = data.df$prof
    bugsData$numProfiles = 102
    bProfileEffect = TRUE

    if(genInits){
      bugsInits$beta.mu = 0
      bugsInits$beta.tau = 1
      #bugsInits$beta.locus = rep(0, bugsData$numProfiles - 1)
    }
  }else{
    bProfileEffect = FALSE
  }

  if(!is.na(m["dye"])){
    bugsData$dye = data.df$dye
    bugsData$numDyes = 4
    bDyeEffect = TRUE

    if(genInits){
      bugsInits$gamma.mu = 0
      bugsInits$gamma.tau = 1
      #bugsInits$gamma.locus = rep(0, bugsData$numDyes - 1)
    }

  }else{
    bDyeEffect = FALSE
  }

  if(!is.na(m["X"])){
    bugsData$X = data.df$X
    bDoseEffect = TRUE
  }else{
    bDoseEffect = FALSE
  }

  if(!is.na(m["V"])){
    bugsData$profile = data.df$prof
    bVarEffect = TRUE
  }else{
    bVarEffect = FALSE
  }


  return(
    list(
      bugsData = bugsData,
      bugsInits = bugsInits,
      bugsModelString = buildModel(responseDist = responseDist,
                                   bLocusEffect = bLocusEffect,
                                   bProfileEffect = bProfileEffect,
                                   bDyeEffect = bDyeEffect,
                                   bDoseEffect = bDoseEffect,
                                   bVarEffect = bVarEffect),
      responseDist = responseDist,
      modelFormula = form,
      effects = list(bLocusEffect = bLocusEffect,
                     bProfileEffect = bProfileEffect,
                     bDyeEffect = bDyeEffect,
                     bDoseEffect = bDoseEffect,
                     bVarEffect = bVarEffect)
    )
  )
}
