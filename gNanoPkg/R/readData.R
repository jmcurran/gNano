#' readData
#'
#' rearrange the BUGSdata into a data.frame
#'
#' @return a data.frame
#' @export
#'
#' @examples
#' data.df = readData()
readData = function(){
  bugsData = readRawData()
  obs = NULL
  loc = NULL
  prof = NULL
  dye = NULL
  X = NULL

  for(allele in 1:2){
    for(locus in 1:bugsData$numLoci){
      loc = c(loc, rep(locus, 102))
      obs = c(obs, bugsData$P[,locus,allele])
      prof = c(prof, 1:102)
      dye = c(dye, bugsData$profileDyes[,locus,allele])
      X = c(X, log(bugsData$X[,locus]))
    }
  }
  freq.df = data.frame(loc = as.factor(loc[!is.na(obs)]),
                       prof = as.factor(prof[!is.na(obs)]),
                       obs = exp(obs[!is.na(obs)]),
                       dye = as.factor(dye[!is.na(obs)]),
                       X = X[!is.na(obs)])

  return(freq.df)
}
