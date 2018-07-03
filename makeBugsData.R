library(tidyverse)

makeBUGSdata = function(LocusNames = c("K1", "M1", "R1", "K2", "Y1", "M2", "R2", "Y2", "Y3", "R3", "R4", "Y4", 
                                       "Y5", "R5", "Y6", "H1", "B1", "S1", "V1", "R6", "Y7", "Y8", "R7", "S2", 
                                       "R8", "Y9", "S3", "R9", "Y10", "R10", "Y11")){
  #reads in file
  rawData = read.csv(here("Duncan/CH_GNano_ForR.csv"), stringsAsFactors = FALSE)
  
  #array of locus names
  allLoci = c("K1", "M1", "R1", "K2", "Y1", "M2", "R2", "Y2", "Y3", "R3", "R4", "Y4", 
              "Y5", "R5", "Y6", "H1", "B1", "S1", "V1", "R6", "Y7", "Y8", "R7", "S2", 
              "R8", "Y9", "S3", "R9", "Y10", "R10", "Y11")
  
  Loci = match.arg(LocusNames, allLoci, several.ok = TRUE)
  numLoci = length(Loci)
  
  ## Filter the rawData down to just the loci we want to use
  rawData = rawData %>% 
    filter(LocusName %in% Loci)
  
  #Number of samples
  numSamples =  max(rawData[ ,2])
  
  
  
  usedDyes = rawData %>% 
    select(Dye1No,Dye2No) %>% 
    distinct() %>% 
    gather() %>% 
    pull(value) %>% 
    unique()
  #gets rid of negative dye values
  usedDyes <- usedDyes[usedDyes > 0]
  
  ## Number of dyes
  numDyes = length(usedDyes)
  
  if(numDyes != 4){
    if(numDyes %in% c(1,3)){
      stop("There are 1 or 3 dyes in use and I am not sure this should happen")
    }
    
    ## This is probably stupidly complicated but I couldn't get it to work
    ## any other way. It writes the code for two mutate operations and 
    ## then evaluates it.
    
    for(i in 1:2){
      mutateStr = paste0("Dye", 
                         i, 
                         "No = recode(Dye", 
                         i, 
                         "No, ", 
                         paste0(sprintf("'%d'  = ", usedDyes), 1:2, collapse = ","), ")") 
      mutateStr = paste0("rawData = rawData %>% mutate(", mutateStr, ")")
      eval(parse(text = mutateStr))
    }
    
  }

  ## Dye names
  dyeNames = c("A", "G", "C", "T")[usedDyes]
  
    
  ## Dye colours
  dyeColours = c("green", "blue", "black", "red")[usedDyes]
  
  #heights formatted into array [sample][locus][allele] 
  profileData = pred = loglik = array(data = NA, dim =c(numSamples, numLoci, 2))
  
  #dyes formatted into array [sample][locus][allele] 
  profileDyes = array(data = NA, dim=c(numSamples, numLoci, 2))
  
  #homozygous formatted into array [sample][locus]
  profileHomozygote = X = alleles_at_locus =  matrix(NA, nrow = numSamples, ncol = numLoci)
  
  
  # library(sqldf)
  library(tidyverse)
  
  loc = 1 ## this counter is to make sure we insert data contiguously
  for(locus in Loci){
     locusData = rawData %>% 
       filter(LocusName == locus)
     
     locusData = locusData %>% 
       mutate(isHom = Allele1 == Allele2,
              allelesAtLocus = ifelse(Allele1 == Allele2, 1, 2),
              dye1 = ifelse(Dye1No %in% 1:4, Dye1No, NA),
              dye2 = ifelse(Dye2No %in% 1:4, Dye2No, NA))
     
     ## set profileHomozygous, alleles_at_Locus, X
     
     sampleNums = locusData %>% pull(SampleNo)
     
     profileHomozygote[sampleNums, loc] = locusData %>% pull(isHom)
     alleles_at_locus[sampleNums, loc]  = locusData %>% pull(allelesAtLocus) 
     X[sampleNums, loc] = ifelse(locusData %>% pull(allelesAtLocus) == 2, 1, 2)
     profileDyes[sampleNums, loc, 1] = locusData %>% pull(dye1)
     
     ## set profileData
     
     profileData[sampleNums, loc, 1] = locusData %>% pull(Height1)
     
     locusData2 = locusData %>% 
       filter(!isHom)
     
     sampleNums2 = locusData2 %>% 
       pull(SampleNo)
     
     profileData[sampleNums2, loc, 2] = locusData2 %>% pull(Height2)
     profileDyes[sampleNums2, loc, 2] = locusData2 %>% pull(dye2)
       
     loc <- loc + 1L
     
     #size of the raw data set
     # rawDataLength = nrow(rawData)
     # 
     # for (r in 1:rawDataLength){
     #   #capturing data from rawData
     #   rawDataProfile = rawData[r, 2]
     #   rawDataLocus = rawData[r, 4]
     #   rawDataDye1 = rawData[r, 7]
     #   rawDataDye2 = rawData[r, 9]
     #   rawDataAlleleHeight1 = rawData[r, 11]
     #   rawDataAlleleHeight2 = rawData[r, 13]
     #   
     #   #determines if homozygous
     #   rawDataAllele1 = rawData[r, 10]
     #   rawDataAllele2 = rawData[r, 12]
     #   
     #   ## This will evaluate as TRUE or FALSE so you don't need the if
     #   profileHomozygote[rawDataProfile, rawDataLocus]  = rawDataAllele1 == rawDataAllele2
     #   X[rawDataProfile, rawDataLocus] =  if(rawDataAllele1 == rawDataAllele2){2}else{1}
     #   if(rawDataAllele1 %in% dyeNames){
     #     alleles_at_locus[rawDataProfile, rawDataLocus] =  if(rawDataAllele1 == rawDataAllele2){1}else{2}
     #   }
     #   
     #   #writing heights to multi-dimensional height array
     #   profileData[rawDataProfile, rawDataLocus, 1] = rawDataAlleleHeight1
     #   if(profileHomozygote[rawDataProfile, rawDataLocus] == FALSE){
     #     profileData[rawDataProfile, rawDataLocus, 2] = rawDataAlleleHeight2
     #   }
     #   
     #   #writing dyes to multi-dimensional height array
     #   profileDyes[rawDataProfile, rawDataLocus, 1] = if(rawDataDye1 %in% 1:4){rawDataDye1}else{NA}
     #   if(!profileHomozygote[rawDataProfile, rawDataLocus]){
     #     profileDyes[rawDataProfile, rawDataLocus, 2] = if(rawDataDye2 %in% 1:4){rawDataDye2}else{NA}
     #   }
     # }
   
     
  }
   
   
  
  
  
  
  locusNo = rawData  %>%  
    filter(LocusName %in% Loci) %>% 
    pull(LocusNo)
  sampleNo = rawData  %>%  
    filter(LocusName %in% Loci) %>% 
    pull(SampleNo)
  
  ## Added by Mikkel
  if (FALSE) {
    str(profileDyes)
    sum(is.na(profileDyes))
    mean(is.na(profileDyes))
    
    all(is.na(profileData) == is.na(profileDyes))
    #dim = c(numSamples, numLoci, 2)
    # remove samples with BOTH NA (one NA is homozygote)
    rm_is <- c() # stupid, but easy and fast enough for this
    for (i in seq_len(dim(profileData)[1])) { # loop over profiles
      #i <- 1
      if (any(apply(profileData[i,Loci,], 1, function(x) all(is.na(x))))) {
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
      locusNo = locusNo[-rm_is]
      sampleNo = sampleNo[-rm_is]
    }
  }
  
  ## Added by Mikkel
  # Per sample, get missing loci:
  #missing_loci <- vector("list", dim(profileData)[1])
  use_loci <- vector("list", dim(profileData)[1])
  for (i in seq_len(dim(profileData)[1])) { # loop over profiles
    #i <- 1
    use_loci[[i]] <- which(!apply(profileData[i,,,drop = FALSE], 1, function(x) all(is.na(x))))
    #use_loci[[i]] <- which(!apply(profileData[i,,], 1, function(x) all(is.na(x))))
  }
  #missing_loci
  
  numSampleLoci = sapply(use_loci, length)
  locOffset = cumsum(numSampleLoci)
  locStart = c(1, locOffset[-length(locOffset)] + 1)
  locEnd = locStart + numSampleLoci - 1
  use_loci = unlist(use_loci)
  
  
  # I want to return every bit of data I KNOW about
  # so all observations and all constants
  bugsData = list(
    numLoci = numLoci,
    numProfiles = numSamples,
    numDyes = numDyes,
    S = 30000,
    profileDyes = profileDyes,
    P = log(profileData),
    use_loci = use_loci,
    locStart = locStart,
    locEnd = locEnd,
    pred = pred, 
    X = X,
    alleles_at_locus = alleles_at_locus,
    loglik = loglik
  )

}
