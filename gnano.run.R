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
  profileData = array(data = NA, dim =c(numSamples, numLoci, 2))
  
  #dyes formatted into array [sample][locus][allele] 
  profileDyes = array(data = NA, dim=c(numSamples, numLoci, 2))
  
  #homozygous formatted into array [sample][locus]
  profileHomozygote = array(data = NA, dim=c(numSamples, numLoci))
  
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
    
    #writing heights to multi-dimensional height array
    profileData[rawDataProfile, rawDataLocus, 1] = rawDataAlleleHeight1
    if(profileHomozygote[rawDataProfile, rawDataLocus] == FALSE){
      profileData[rawDataProfile, rawDataLocus, 2] = rawDataAlleleHeight2
    }
    
    #writing dyes to multi-dimensional height array
    profileDyes[rawDataProfile, rawDataLocus, 1] = rawDataDye1
    if(!profileHomozygote[rawDataProfile, rawDataLocus]){
      profileDyes[rawDataProfile, rawDataLocus, 2] = rawDataDye2
    }
  }
  
  bugsData = list(
    numLoci = numLoci,
    numProfiles = numSamples,
    P = profileData, 
    
    
  )
}