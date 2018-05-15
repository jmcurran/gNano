model{
  #Amplification Priors
  for(locus in 1:numLoci){
    mu.amp[locus] ~ dunif(0.1, 10)
    sigma.sq.amp[locus] ~ dunif(0.01, 1)
    tau.amp[locus] <- 1 / sigma.sq.amp[locus]
  }
  
  #Dye effect priors
  for(f in 1:numDyes){
    mu.dye[f] ~ dunif(0.1, 10)
    sigma.sq.dye[f] ~ dunif(0.01, 1)
    tau.dye[f] <- 1  / sigma.sq.dye[f]
  }
  
  #Lambda prior
  lambda ~ dunif(1, 100)
  
  #Profiles
  for(c in 1:numProfiles){
    #draw a locus amp efficiency for each locus
    for(locus in 1:numLoci){
      A[c,locus] ~ dlnorm(mu.amp[locus], tau.amp[locus])
    }
    
    #draw a dye effect for each fluorophore
    for(f in 1:numDyes){
      D[c, f] ~ dlnorm(mu.dye[f], tau.dye[f])
    }
    
    #choose a template for each contributor
    T[c] ~ dunif(0, S)
    #peak height variance for this profile
    Var[c] <- lambda / T[c]
    Prec[c] <- T[c] / lambda
    
    #generate expected heights for each allele at each locus in each profile
    for(locus in 1:numLoci){ 
      for(allele in 1:alleles_at_locus[c,locus]){
        E[c, locus, allele] = T[c] * A[c, locus] * D[c, profileDyes[c, locus, allele]] * X[locus, allele]
        #model observed heights
        P[c, locus, allele] ~ dlnorm(log(E[c, locus, allele]), Prec[c])
      }
    }
  }
}