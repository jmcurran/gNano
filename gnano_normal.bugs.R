model{
  #Amplification Priors
  Mu.amp ~ dnorm(0, 0.000001)
  for(locus in 1:numLoci){
    mu.amp[locus] ~ dnorm(Mu.amp, 0.000001)
  }
  
  #Dye effect priors
  Mu.dye ~ dnorm(0, 0.000001)
  for(f in 1:numDyes){
    mu.dye[f] ~ dnorm(Mu.dye, 0.000001)
  }
  
  #Lambda prior
  lambda ~ dunif(1, 500)
  
  mu ~ dnorm(0, 0.000001)
  
  
  #Profiles
  for(c in 1:numProfiles){
    #choose a template for each contributor
    T[c] ~ dunif(-log(S), log(S))
    #peak height variance for this profile
    Var[c] <- lambda / exp(T[c])
    Prec[c] <- 1 / Var[c]
    
    
    
    #generate expected heights for each allele at each locus in each profile
    for(locus in use_loci[locStart[c]:locEnd[c]]){ 
      for(allele in 1:alleles_at_locus[c,locus]){
        E[c, locus, allele] = log(X[c, locus]) + T[c] + mu.amp[locus] + mu.dye[profileDyes[c, locus, allele]]  
        #model observed heights
        P[c, locus, allele] ~ dnorm(E[c, locus, allele], Prec[c])T(0, 30000)
        pred[c, locus, allele] ~ dnorm(E[c, locus, allele], Prec[c])T(0, 30000)
        loglik[c, locus, allele] <- log(dnorm(P[c, locus, allele], E[c, locus, allele], Prec[c]))
      }
    }
  }
}