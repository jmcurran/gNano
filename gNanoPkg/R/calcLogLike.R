#' Calculate the log-likelihood
#'
#' @param results.df
#' @param responseDist
#'
#' @return
#' @export
#'
#' @examples
calcLogLik = function(results.df, responseDist = c("g", "ln", "sn"), n.iter = 1000){

  responseDist = match.arg(responseDist)

  gNano.df = readData()
  nObs = length(gNano.df$obs)

  if(responseDist == "g"){
    sr.df = results.df %>%
      select(c(starts_with("shape"), starts_with("rate")))

    sr.df = sr.df %>%
      pivot_longer(cols = everything(), names_to = "parameter")



    ## test to see if shape and rate have per obs. values
    ## i.e do we have shape[1], shape[2], ... , shape[4069]
    ## this will be true for all models above g-0

    bBrackets = (nrow(sr.df) / (2 * n.iter)) > 1

    if(bBrackets){
      sr.df = sr.df %>%
        mutate(obs = as.numeric(gsub("^(shape|rate)([.]|\\[)([0-9]+)([.]|\\])$", "\\3", parameter)))

      sr.df = sr.df %>%
        mutate(parameter = gsub("^(shape|rate)([.]|\\[)([0-9]+)([.]|\\])$", "\\1", parameter))
    }

    shape.df = sr.df %>%
      filter(parameter == "shape")

    rate.df = sr.df %>%
      filter(parameter == "rate")

    sr.df = if(nrow(shape.df) == n.iter){
      data.frame(
        obs = rep(1:nObs, n.iter),
        y = rep(gNano.df$obs, n.iter),
        shape = shape.df$value,
        rate = rate.df$value
      )
    }else{
      data.frame(
        obs = shape.df$obs,
        y = rep(gNano.df$obs, n.iter),
        shape = shape.df$value,
        rate = rate.df$value
      )
    }


    sr.df = sr.df %>%
      mutate(rep = rep(1:n.iter, rep(nObs, n.iter)))

    ll.df = sr.df %>%
      group_by(rep) %>%
      summarise(ll = sum(dgamma(y, shape, rate, log = TRUE)))

  }else if(responseDist == "ln"){
    ms.df = results.df %>%
      select(matches("^(mu|tau)(([.]|\\[)[0-9]+([.]|\\]))?$"))

    ms.df = ms.df %>%
      pivot_longer(cols = everything(), names_to = "parameter") %>%
      mutate(parameter = tolower(parameter))

    nObs = length(gNano.df$obs)

    bBrackets = any(grepl("\\[", ms.df$parameter))

    if(bBrackets){
      ms.df = ms.df %>%
        mutate(obs = as.numeric(gsub("^(mu|tau)([.]|\\[)([0-9]+)([.]|\\])$", "\\3", parameter)))

      ms.df = ms.df %>%
        mutate(parameter = gsub("^(mu|tau)([.]|\\[)([0-9]+)([.]|\\])$", "\\1", parameter))
    }

    mu.df = ms.df %>%
      filter(parameter == "mu")

    sigma.df = ms.df %>%
      filter(parameter == "tau") %>%
      mutate(sigma = 1 / sqrt(value)) %>%
      select(-value)

    ms.df = if(nrow(mu.df) == n.iter && nrow(sigma.df) == n.iter){
      data.frame(
        obs = rep(1:nObs, n.iter),
        y = rep(gNano.df$obs, n.iter),
        mu = rep(mu.df$value, rep(nObs, n.iter)),
        sigma = rep(sigma.df$sigma, rep(nObs, n.iter))
      )
    }else if(nrow(mu.df) > n.iter && nrow(sigma.df) == n.iter){
      data.frame(
        obs = rep(1:nObs, n.iter),
        y = rep(gNano.df$obs, n.iter),
        mu = mu.df$value,
        sigma = rep(sigma.df$sigma, rep(nObs, n.iter))
      )
    }else if(nrow(mu.df) == n.iter && nrow(sigma.df) > n.iter){
      data.frame(
        obs = rep(1:nObs, n.iter),
        y = rep(gNano.df$obs, n.iter),
        mu = rep(mu.df$value, rep(nObs, n.iter)),
        sigma = sigma.df$sigma
      )
    }else{
      data.frame(
        obs = mu.df$obs,
        y = rep(gNano.df$obs, n.iter),
        mu = mu.df$value,
        sigma = sigma.df$sigma
      )
    }

    ms.df = ms.df %>%
      mutate(rep = rep(1:n.iter, rep(nObs, n.iter)))

    ll.df = ms.df %>%
      group_by(rep) %>%
      summarise(ll = sum(dlnorm(y, mu, sigma, log = TRUE)))
  }else{ ## "sn"
    lss.df = results.df %>%
      pivot_longer(cols = everything(), names_to = "parameter")

     lss.df = lss.df %>%
      mutate(parameter = gsub("^(location|skew|scale)(([.]|\\[)([0-9]+)([.]|\\]))?$", "\\1", parameter))

    location.df = lss.df %>%
      filter(parameter == "location")

    location.df = lss.df[lss.df$parameter == "location",]
    scale.df = lss.df[lss.df$parameter == "scale",]
    skew.df = lss.df[lss.df$parameter == "skew",]

    lss.df = if(nrow(location.df) == n.iter && nrow(scale.df) == n.iter){
      lss.df = data.frame(
        obs = rep(1:nObs, n.iter),
        y = rep(log(gNano.df$obs), n.iter),
        location = rep(location.df$value, rep(nObs, n.iter)),
        scale = rep(scale.df$value, rep(nObs, n.iter)),
        skew = rep(skew.df$value, rep(nObs, n.iter))
      )
    }else if(nrow(location.df) > n.iter && nrow(scale.df) == n.iter){
      lss.df = data.frame(
        obs = rep(1:nObs, n.iter),
        y = rep(log(gNano.df$obs), n.iter),
        location = location.df$value,
        scale = rep(scale.df$value, rep(nObs, n.iter)),
        skew = rep(skew.df$value, rep(nObs, n.iter))
      )

    }else if(nrow(location.df) == n.iter && nrow(scale.df) > n.iter){
      lss.df = data.frame(
        obs = rep(1:nObs, n.iter),
        y = rep(log(gNano.df$obs), n.iter),
        location = rep(location.df$valuerep(nObs, n.iter)),
        scale = scale.df$value,
        skew = rep(skew.df$value, rep(nObs, n.iter))
      )
    }else{
      lss.df = data.frame(
        obs = rep(1:nObs, n.iter),
        y = rep(log(gNano.df$obs), n.iter),
        location = location.df$value,
        scale = scale.df$value,
        skew = rep(skew.df$value, rep(nObs, n.iter))
      )
    }

    lss.df = lss.df %>%
      mutate(rep = rep(1:n.iter, rep(nObs, n.iter)))

    ll.df = lss.df %>%
      group_by(rep) %>%
      summarise(ll = sum(.dsn(y, xi = location, omega = scale, alpha = skew, log = TRUE)))
  }

  return(ll.df)
}
