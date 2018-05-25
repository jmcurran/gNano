model{
  #Lambda prior
  lambda ~ dunif(1, 500)
  
  #Profiles
  for(c in 1:numProfiles){    
    #choose a template for each contributor
    T[c] ~ dunif(0, S)
    #peak height variance for this profile
    Var[c] <- lambda / T[c]
    Prec[c] <- T[c] / lambda
   
    #generate expected heights for each allele at each locus in each profile
    for(locus in use_loci[locStart[c]:locEnd[c]]){ 
      for(allele in 1:alleles_at_locus[c,locus]){
        E[c, locus, allele] = T[c] * X[locus, allele]
        #model observed heights
        #LNmu <- log(E[c, locus, allele])+Var[c]
        P[c, locus, allele] ~ dlnorm(log(E[c, locus, allele])+Var[c], Prec[c])T(0, 30000)
        pred[c, locus, allele] ~ dlnorm(log(E[c, locus, allele])+Var[c], Prec[c])T(0, 30000)
        loglik[c, locus, allele] <- log(dlnorm(P[c, locus, allele], log(E[c, locus, allele])+Var[c], Prec[c]))
      }
    }
  }
}