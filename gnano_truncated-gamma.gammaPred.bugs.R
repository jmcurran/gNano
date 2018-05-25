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
  
  #Dye effect priors
  for(f in 1:numDyes){
    tau.dye[f] ~ dgamma(0.01, 0.01)
    sigma.sq.dye[f] <- 1 / tau.dye[f]
  }
  for(f in 1:(numDyes-1)){
    mu.dye[f] ~ dunif(0, 10)
  }
  #final dye effect is calculated so that geomean = 1
  mu.dye[numDyes] <- 1/prod(mu.dye[1:(numDyes-1)])
  
  #Lambda prior
  lambda ~ dunif(1, 500)
  
  #Profiles
  for(c in 1:numProfiles){
    #draw a locus amp efficiency for each locus
    for(locus in 1:numLoci){
      A[c,locus] ~ dlnorm(log(mu.amp[locus]), tau.amp[locus])
    }
    
    #draw a dye effect for each fluorophore
    for(f in 1:numDyes){
      D[c, f] ~ dlnorm(log(mu.dye[f]), tau.dye[f])
    }
    
    #choose a template for each contributor
    T[c] ~ dunif(0, S)

    #generate expected heights for each allele at each locus in each profile
    for(locus in use_loci[locStart[c]:locEnd[c]]){ 
      for(allele in 1:alleles_at_locus[c,locus]){
        E[c, locus, allele] = T[c] * A[c, locus] * D[c, profileDyes[c, locus, allele]] * X[locus, allele]
        #model observed heights
        GammaRate[c, locus, allele] <- E[c, locus, allele]*(exp((3*lambda)/(2*T[c]))-1)
        GammaShape[c, locus, allele] <- E[c, locus, allele]/GammaRate[c, locus, allele]+1
        P[c, locus, allele] ~ dgamma(GammaShape[c, locus, allele], GammaRate[c, locus, allele])T(0, 30000)
        pred[c, locus, allele] ~ dgamma(GammaShape[c, locus, allele], GammaRate[c, locus, allele])T(0, 30000)
        loglik[c, locus, allele] <- log(dgamma(P[c, locus, allele], GammaShape[c, locus, allele], GammaRate[c, locus, allele]))
      }
    }
  }
}