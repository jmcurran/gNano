model{
  #Amplification Priors
  for(locus in 1:numLoci){
    tau.amp[locus] ~ dgamma(0.01, 0.01)
    sigma.sq.amp[locus] <- 1 / tau.amp[locus]
    mu.amp[locus] ~ dnorm(0, 0.000001)
  }
  
  #Dye effect priors
  for(f in 1:numDyes){
    tau.dye[f] ~ dgamma(0.01, 0.01)
    sigma.sq.dye[f] <- 1 / tau.dye[f]
    mu.dye[f] ~ dnorm(0, 0.000001)
  }
  
  #final dye effect is calculated so that geomean = 1
  # mu.dye[numDyes] <- 1/prod(mu.dye[1:(numDyes-1)])
  
  #Lambda prior
  lambda ~ dunif(1, 500)
  
  #Profiles
  for(c in 1:numProfiles){
    #draw a locus amp efficiency for each locus
    for(locus in 1:numLoci){
      A[c,locus] ~ dnorm(mu.amp[locus], tau.amp[locus])
    }
    
    #draw a dye effect for each fluorophore
    for(f in 1:numDyes){
      D[c, f] ~ dnorm(mu.dye[f], tau.dye[f])
    }
    
    #choose a template for each contributor
    T[c] ~ dunif(0, S)
    #peak height variance for this profile
    Var[c] <- lambda / T[c]
    Prec[c] <- T[c] / lambda
    
    
    #generate expected heights for each allele at each locus in each profile
    for(locus in use_loci[locStart[c]:locEnd[c]]){ 
      for(allele in 1:alleles_at_locus[c,locus]){
        log(E[c, locus, allele]) = X[locus, allele]*T[c] + A[c, locus] + D[c, profileDyes[c, locus, allele]]  
        #model observed heights
        #LNmu <- log(E[c, locus, allele])+Var[c]
        P[c, locus, allele] ~ dlnorm(E[c, locus, allele], Prec[c])T(0, 30000)
        #pred[c, locus, allele] ~ dnorm(E[c, locus, allele], Prec[c])T(0, 30000)
        #loglik[c, locus, allele] <- log(dnorm(P[c, locus, allele], E[c, locus, allele], Prec[c]))
      }
    }
  }
}