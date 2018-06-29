calculateWAIC <- function(sim.sample, saveDir, bugsData) {
  #library required for WAIC (also requires matrixstats)
  library(loo)
  NumberLoci <- bugsData$numLoci
  profileData <- bugsData$P
  profileDyes <- bugsData$profileDyes
  NumberSamples <- dim(profileData)[[1]]
  alleles_at_locus <- bugsData$alleles_at_locus
  Number_iterations <- length(sim.sample[[1]][, "loglik[1,1,1]"])
  number_chains <- length(sim.sample[, "loglik[1,1,1]"])
  
  #capture the log-likleihood data in an array
  logLik.unif <-
    array(NA,
          dim = c(
            number_chains,
            Number_iterations,
            NumberSamples,
            NumberLoci,
            2
          ))
  for (chain in 1:number_chains) {
    for (sample in 1:NumberSamples) {
      for (locus in 1:NumberLoci) {
        for (allele in 1:alleles_at_locus[sample, locus]) {
          if (!is.na(profileData[sample, locus, allele])) { 
            logLik.unif[chain, , sample, locus, allele] <-
              sim.sample[[chain]][, paste("loglik[", sample, ",", locus, ",", allele, "]", sep ="")]
          }    
        }
      }
    }
  }
  
  #calculates WAIC, again summing log-likelihoods across alleles at a locus
  waic.unif <- waic(apply(logLik.unif, c(1, 2, 3), sum, na.rm = TRUE))
  
  #comaprison of waic values
  print("WAIC for uniform model")
  waic.unif
  
  #saves to text file
  pathToFile <- paste(saveDir, "WAICresults.txt", sep = "")
  write.table(waic.unif[[1]],
              file = pathToFile,
              col.names = TRUE,
              quote = F)
}
