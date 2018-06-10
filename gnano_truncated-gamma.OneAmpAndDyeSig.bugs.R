model{
  #Amplification Priors
  tau.amp ~ dgamma(0.01, 0.01)
  sigma.sq.amp <- 1 / tau.amp
  for(locus in 1:(numLoci-1)){
    mu.amp[locus] ~ dunif(0, 10)
  }
  #final amp effect is calculated so that geomean = 1
  mu.amp[numLoci] <- 1/prod(mu.amp[1:(numLoci-1)])
  
  #Dye effect priors
  tau.dye ~ dgamma(0.01, 0.01)
  sigma.sq.dye <- 1 / tau.dye
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
      A[c,locus] ~ dlnorm(log(mu.amp[locus]), tau.amp)
    }
    
    #draw a dye effect for each fluorophore
    for(f in 1:numDyes){
      D[c, f] ~ dlnorm(log(mu.dye[f]), tau.dye)
    }
    
    #choose a template for each contributor
    T[c] ~ dunif(0, S)
    #peak height variance for this profile
    Var[c] <- lambda / T[c]
    Prec[c] <- T[c] / lambda
    
    
    #generate expected heights for each allele at each locus in each profile
    for(locus in use_loci[locStart[c]:locEnd[c]]){ 
      for(allele in 1:alleles_at_locus[c,locus]){
        E[c, locus, allele] = T[c] * A[c, locus] * D[c, profileDyes[c, locus, allele]] * X[locus, allele]
        #model observed heights
        #LNmu <- log(E[c, locus, allele])+Var[c]
        P[c, locus, allele] ~ dlnorm(log(E[c, locus, allele])+Var[c], Prec[c])T(0, 30000)
        pred[c, locus, allele] ~ dlnorm(log(E[c, locus, allele])+Var[c], Prec[c])T(0, 30000)
        loglik[c, locus, allele] <- log(dlnorm(P[c, locus, allele], log(E[c, locus, allele])+Var[c], Prec[c]))
      }
    }
  }
}