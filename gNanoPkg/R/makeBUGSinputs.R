makeBUGSdata = function(form, data.df, responseDist = c("gamma", "normal")){
  if(!is_formula(form)){
    stop("form must be a formula")
  }

  responseDist = match.arg(reponseDist)

  vars = attr(terms(form), "term.labels")

  if(!all(grepl("locus|dye|profile|X", vars))){
    mismatch = vars[-grep("locus|dye|profile|X", vars)]
    paste0("The variables : ", mismatch, " don't match the list of allowable variables" )
  }

  ## lazy - always assume y and number of observations are in the reponse
  ## will need to hand logs though
  bugsData = list(y = ifelse(responseDist == "gamma", data.df$obs, log(data.df,obs)),
                  N = length(data.df$obs))

  m = match(c("locus", "profile", "dye", "X"), vars)
  names(m) = c("locus", "profile", "dye", "X")

  if(!is.na(m["locus"])){
    bugsData$locus = data.df$loc
    bugsData$numLoci = 31 ## this may not be true
  }

  if(!is.na(m["profile"])){
    bugsData$profile = data.df$profile
    bugsData$numProfiles = 102
  }

  if(!is.na(m["dye"])){
    bugsData$dye = data.df$dye
    bugsData$numDyes = 4
  }

  if(!is.na(m["X"])){
    bugsData$X = data.df$X
  }

  return(bugsData)
}
