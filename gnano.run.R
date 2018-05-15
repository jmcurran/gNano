library(here)

makeBUGSdata = function(){
  #reads in file
  rawData = read.csv(here("Duncan/CH_GNano_ForR.csv"), stringsAsFactors = FALSE)
  
  #array of locus names
  LocusNames = c("K1", "M1", "R1", "K2", "Y1", "M2", "R2", "Y2", "Y3", "R3", "R4", "Y4", 
                 "Y5", "R5", "Y6", "H1", "B1", "S1", "V1", "R6", "Y7", "Y8", "R7", "S2", 
                 "R8", "Y9", "S3", "R9", "Y10", "R10", "Y11")
  numLoci = length(LocusNames)
  
  #Number of samples
  numSamples =  max(rawData[ ,2])
  
  #Dye names
  dyeNames = c("A", "G", "C", "T")
  
  #number of dyes
  numDyes = length(dyeNames)
  
  #dye colours
  dyeColours = c("green", "blue", "black", "red")
  
  #heights formatted into array [sample][locus][allele] 
  profileData = pred = array(data = NA, dim =c(numSamples, numLoci, 2))
  
  #dyes formatted into array [sample][locus][allele] 
  profileDyes = array(data = NA, dim=c(numSamples, numLoci, 2))
  
  #homozygous formatted into array [sample][locus]
  profileHomozygote = X = alleles_at_locus =  matrix(NA, nrow = numSamples, ncol = numLoci)
  
  
  # library(sqldf)
  # library(tidyverse)
  # 
  # sampleNames = rawData %>% 
  #   select(Sample) %>% 
  #   unique() #%>% 
  # 
  
  
  #size of the raw data set
  rawDataLength = nrow(rawData)
  
  for (r in 1:rawDataLength){
    #capturing data from rawData
    rawDataProfile = rawData[r, 2]
    rawDataLocus = rawData[r, 4]
    rawDataDye1 = rawData[r, 7]
    rawDataDye2 = rawData[r, 9]
    rawDataAlleleHeight1 = rawData[r, 11]
    rawDataAlleleHeight2 = rawData[r, 13]
    
    #determines if homozygous
    rawDataAllele1 = rawData[r, 10]
    rawDataAllele2 = rawData[r, 12]
    
    ## This will evaluate as TRUE or FALSE so you don't need the if
    profileHomozygote[rawDataProfile, rawDataLocus]  = rawDataAllele1 == rawDataAllele2
    X[rawDataProfile, rawDataLocus] =  if(rawDataAllele1 == rawDataAllele2){2}else{1}
    if(rawDataAllele1 %in% dyeNames){
      alleles_at_locus[rawDataProfile, rawDataLocus] =  if(rawDataAllele1 == rawDataAllele2){1}else{2}
    }
    
    #writing heights to multi-dimensional height array
    profileData[rawDataProfile, rawDataLocus, 1] = rawDataAlleleHeight1
    if(profileHomozygote[rawDataProfile, rawDataLocus] == FALSE){
      profileData[rawDataProfile, rawDataLocus, 2] = rawDataAlleleHeight2
    }
    
    #writing dyes to multi-dimensional height array
    profileDyes[rawDataProfile, rawDataLocus, 1] = if(rawDataDye1 %in% 1:4){rawDataDye1}else{NA}
    if(!profileHomozygote[rawDataProfile, rawDataLocus]){
      profileDyes[rawDataProfile, rawDataLocus, 2] = if(rawDataDye2 %in% 1:4){rawDataDye2}else{NA}
    }
  }
  
  ## Added by Mikkel
  if (TRUE) {
    str(profileDyes)
    sum(is.na(profileDyes))
    mean(is.na(profileDyes))
    
    all(is.na(profileData) == is.na(profileDyes))
    #dim = c(numSamples, numLoci, 2)
    # remove samples with BOTH NA (one NA is homozygote)
    rm_is <- c() # stupid, but easy and fast enough for this
    for (i in seq_len(dim(profileData)[1])) { # loop over profiles
      #i <- 1
      if (any(apply(profileData[i,,], 1, function(x) all(is.na(x))))) {
        rm_is <- c(rm_is, i)
      }
    }
    rm_is
    length(rm_is)
    
    if (length(rm_is) > 0L) {
      profileData <- profileData[-rm_is,,]
      profileDyes <- profileDyes[-rm_is,,]
      X <- X[-rm_is, ]
      numSamples <- numSamples - length(rm_is)
      alleles_at_locus <- alleles_at_locus[-rm_is, ]
    }
  }
  
  # I want to return every bit of data I KNOW about
  # so all observations and all constants
  bugsData = list(
    numLoci = numLoci,
    numProfiles = numSamples,
    numDyes = numDyes,
    S = 30000,
    profileDyes = profileDyes,
    P = profileData, 
    pred = pred, 
    X = X,
    alleles_at_locus = alleles_at_locus
  )
}

library(rjags)
bugsData = makeBUGSdata()
bugsFile = here("gnano.bugs.R")

## compile the model
sim = jags.model(file = bugsFile, data = bugsData)

## do a bit of burn in - no idea what is sufficient at this point
system.time(update(sim, 100000))

## What parameters are we interested?

parameters = c("pred")
sim.sample = coda.samples(model = sim, variable.names = parameters,
                          n.iter = 50000, thin = 50)
stats = summary(sim.sample)
#stats
#plot(sim.sample)
