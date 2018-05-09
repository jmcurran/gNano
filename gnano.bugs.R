model{
  #Amplification Priors
  for(locus in 1:numLoci){
    mu.amp[locus] ~ dunif(0.1, 10)
    sigma.sq.amp[locus] ~ dunif(0.01, 1)
    tau.amp[locus] <- 1 / sigma.sq.amp[locus]
  }
  
  #Dye effect priors
  for(f in 1:numFlurophores){
    mu.dye[f] ~ dunif(0.1, 10)
    sigma.sq.dye[f] ~ dunif(0.01, 1)
  }
  
  #Lambda prior
  lambda ~ dunif(1, 100)
  
  #Profiles
  for(c in 1:numProf){
    #draw a locus amp efficiency for each locus
    for(locus in 1:numLoci){
      A[c,locus] ~ dlnorm(mu.amp[locus], tau.amp[locus])
    }
    
    #draw a dye effect for each fluorophore
    for(f in 1:numFluorophores){
      D[c, f] ~ dlnorm(mu.dye[locus], tau.dye[locus])
    }
    
    #choose a template for each contributor
    T[c] ~ dunif(0, S)
    #peak height variance for this profile
    Var[c] <- lambda / T[c]
    
    #generate expected heights for each allele at each locus in each profile
    for(locus in 1:numLoci){ 
      for(allele in 1:2){
        E[c, locus, allele] = T[c] * A[c, locus] * D[c, allele_to_fluorophore[Allele[c,locus, allele]]] * X[c,locus,allele]
        #model observed heights
        P[c, locus, allele] ~ dlnorm(ln(E[c, l, a]), sqrt(Var[c]))
      }
    }
  }
}