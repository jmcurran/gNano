model{
  #Amplification Priors
  for(locus in 1:numLoci){
    tau.amp[locus] ~ dgamma(0.01, 0.01)
    sigma.sq.amp[locus] <- 1 / tau.amp[locus]
  }
  for(locus in 1:(numLoci-1)){
    mu.amp[locus] ~ dunif(0, 10)
  }
  #final amp effect is calculated so that geomean = 1
  mu.amp[numLoci] <- 1/prod(mu.amp[1:(numLoci-1)])

  #Lambda prior
  lambda ~ dunif(1, 500)
  
  #Profiles
  for(c in 1:numProfiles){
    #draw a locus amp efficiency for each locus
    for(locus in 1:numLoci){
      A[c,locus] ~ dlnorm(log(mu.amp[locus]), tau.amp[locus])
    }
    
    #choose a template for each contributor
    T[c] ~ dunif(0, S)
    #peak height variance for this profile
    Var[c] <- lambda / T[c]
    Prec[c] <- T[c] / lambda
   
    
    #generate expected heights for each allele at each locus in each profile
    for(locus in use_loci[locStart[c]:locEnd[c]]){ 
      for(allele in 1:alleles_at_locus[c,locus]){
        E[c, locus, allele] = T[c] * A[c, locus] * X[locus, allele]
        #model observed heights
        P[c, locus, allele] ~ dlnorm(log(E[c, locus, allele]), Prec[c])T(0, 30000)
        pred[c, locus, allele] ~ dlnorm(log(E[c, locus, allele]), Prec[c])T(0, 30000)
      }
    }
  }
}