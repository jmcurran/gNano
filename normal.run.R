#load required libraries
library(here)
library(rjags)

#loads the make bugs data method
source("makeBugsData.R")

#creates datafile
bugsData = makeBUGSdata(c("R1", "K2"))
nChains = 1

bugsFile = here("gnano_normal.bugs.R")

## compile the model
system.time({sim = jags.model(file = bugsFile,
                 data = bugsData,
                 n.chains = nChains)})

## do a bit of burn in - no idea what is sufficient at this point
system.time(update(sim, 100000))
#system.time(update(sim, 100))

## The parameters are we interested
parameters = c("alpha.amp",
               "beta.dye",
               "pred",
               "lambda",
               "mu.dye",
               "mu.amp",
               "Mu.dye",
               "Mu.amp",
               "T",
               "loglik")

## run the model
system.time({
sim.sample = coda.samples(
  model = sim,
  variable.names = parameters,
  n.iter = 50000,
  thin = 50
)})

#creates graphs
#source("gNano.graphs.R")
#saveDir <- here("graphs - uniform//")
#saveDir <- here(paste("graphs - ", model_descriptor, "//", sep = ""))
#createGraphs(sim.sample, saveDir, bugsData)

#run the WAIC
#source("gNano.WAIC.R")
#calculateWAIC(sim.sample, saveDir, bugsData)

## save sample
#save(sim.sample, file = paste0(saveDir, "sim.sample.Rda"))

#save summary stats
simSummary <- summary(sim.sample)
#pathToFile <- paste(saveDir, "simStatsSummary.txt", sep = "")
#write.table(simSummary[[1]],
#             file = pathToFile,
#             col.names = TRUE,
#             quote = F)
# pathToFile <- paste(saveDir, "simStatsSummaryQuantiles.txt", sep = "")
# write.table(simSummary[[2]],
#             file = pathToFile,
#             col.names = TRUE,
#             quote = F)

i = grep("pred" ,rownames(simSummary$statistics))
yhat = simSummary$statistics[i,1]
obs = c(bugsData$P[,1,1], bugsData$P[,2,1], bugsData$P[,1,2], bugsData$P[,2,2])

obs = NULL
for(allele in 1:2)
  for(locus in 1:bugsData$numLoci){
    obs = c(obs, bugsData$P[,locus,allele])
  }

obs = obs[!is.na(obs)]

plot(obs,yhat)

i = grep("mu.dye", colnames(sim.sample[[1]]))
plot(sim.sample[[1]][,i])


################# 
#DUNCAN - this is where we were when you left 4/7/18
#

bugsData = makeBUGSdata()

obs = NULL
loc = NULL
prof = NULL
dye = NULL
X = NULL
for(allele in 1:2){
  for(locus in 1:bugsData$numLoci){
    loc = c(loc, rep(locus, 102))
    obs = c(obs, bugsData$P[,locus,allele])
    prof = c(prof, 1:102)
    dye = c(dye, bugsData$profileDyes[,locus,allele])
    X = c(X, log(bugsData$X[,locus]))
  }
}
freq.df = data.frame(loc = as.factor(loc[!is.na(obs)]),
                     prof = as.factor(prof[!is.na(obs)]),
                     obs = obs[!is.na(obs)],
                     dye = as.factor(dye[!is.na(obs)]),
                     X = X[!is.na(obs)])

## This model has a locus effect, and it treats the profile as a random effect
fit = aov(obs ~ loc + dye + Error(prof) + X, data = freq.df)
summary(fit) ## and this shows we're well-justified to say there is a locus effect.
rm(obs)
rm(loc)
rm(prof)
rm(dye)
rm(X)

#################### the simple bugs version with only locus and profile effects #######################

bugsData = list(y = freq.df$obs, locus = freq.df$loc, profile = freq.df$prof, N = length(freq.df$obs),
                numLoci = 31, numProfiles = 102)
nChains = 1

bugsFile = here("simple.bugs.R")

## compile the model
system.time({sim = jags.model(file = bugsFile,
                              data = bugsData,
                              n.chains = nChains)})
update(sim, 10000)
parameters = c("Mu", "alpha.locus", "mu")
#parameters = c("Mu","sigma")
sim.sample = coda.samples(sim, parameters, n.iter = 1000)
simSummary = summary(sim.sample)

i = grep("^(alpha).*$", rownames((simSummary$statistics)))

library(Hmisc)
errbar(1:31, simSummary$quantiles[i,3], simSummary$quantiles[i,1], simSummary$quantiles[i,5],
       xlab = "Locus", ylab = "Locus Effects")

i = grep("^(mu).*$", rownames((simSummary$statistics)))
fitted = simSummary$statistics[i,1] ## means
plot(fitted~bugsData$y, xlab = "Observed", ylab = "Fitted")
abline(c(0,1), col = "red")

################# the simple bugs version with locus, profile and dye effects ##########################

bugsData = list(y = freq.df$obs, locus = freq.df$loc, dye = freq.df$dye, profile = freq.df$prof, N = length(freq.df$obs),
                numLoci = 31, numProfiles = 102, numDyes = 4)
nChains = 1

bugsFile = here("simple.bugs.withDye.R")

## compile the model
system.time({sim = jags.model(file = bugsFile,
                              data = bugsData,
                              n.chains = nChains)})
update(sim, 10000)
parameters = c("Mu", "alpha.locus", "gamma.dye", "mu")
sim.sample = coda.samples(sim, parameters, n.iter = 1000)
simSummary = summary(sim.sample)

#graphing the locus effects
i = grep("^(alpha).*$", rownames((simSummary$statistics)))
library(Hmisc)
errbar(1:31, simSummary$quantiles[i,3], simSummary$quantiles[i,1], simSummary$quantiles[i,5],
       xlab = "Locus", ylab = "Locus Effects")

#graphing the dye effects
i = grep("^(gamma).*$", rownames((simSummary$statistics)))
library(Hmisc)
numDyes = 4
errbar(1:numDyes, simSummary$quantiles[i,3], simSummary$quantiles[i,1], simSummary$quantiles[i,5],
       xlab = "Dye", ylab = "Dye Effects",
       axes = FALSE)
axis(2)
axis(1, at = 1:4)
box()

#graphing obs vs expected
i = grep("^(mu).*$", rownames((simSummary$statistics)))
fitted = simSummary$statistics[i,1] ## means
plot(fitted~bugsData$y, xlab = "Observed", ylab = "Fitted")
abline(c(0,1), col = "red")

         
################# the medium bugs version with locus, profile and dye effects ##########################

bugsData = list(y = freq.df$obs, locus = freq.df$loc, dye = freq.df$dye, profile = freq.df$prof, N = length(freq.df$obs),
                numLoci = 31, numProfiles = 102, numDyes = 4)
nChains = 1

bugsFile = here("medium.bugs.withDye.R")

## compile the model
system.time({sim = jags.model(file = bugsFile,
                              data = bugsData,
                              n.chains = nChains)})
update(sim, 10000)
parameters = c("Mu", "alpha.mu", "alpha.sigma", "alpha.locus", "gamma.dye", "mu")
sim.sample = coda.samples(sim, parameters, n.iter = 1000)
simSummary = summary(sim.sample)

#graphing the locus effects
i = grep("^(alpha\\.locus).*$", rownames((simSummary$statistics)))
library(Hmisc)
errbar(1:31, simSummary$quantiles[i,3], simSummary$quantiles[i,1], simSummary$quantiles[i,5],
       xlab = "Locus", ylab = "Locus Effects")

#graphing the dye effects
i = grep("^(gamma).*$", rownames((simSummary$statistics)))
library(Hmisc)
numDyes = 4
errbar(1:numDyes, simSummary$quantiles[i,3], simSummary$quantiles[i,1], simSummary$quantiles[i,5],
       xlab = "Dye", ylab = "Dye Effects",
       axes = FALSE)
axis(2)
axis(1, at = 1:4)
box()

#graphing obs vs expected
i = grep("^(mu).*$", rownames((simSummary$statistics)))
fitted = simSummary$statistics[i,1] ## means
plot(fitted~bugsData$y, xlab = "Observed", ylab = "Fitted")
abline(c(0,1), col = "red")

results.df = data.frame(fitted = fitted, observed = bugsData$y)

library(ggplot2)
library(tidyverse)
p = results.df %>% ggplot(aes(x = observed, y = fitted)) + geom_point() + geom_abline(slope = 1, intercept = 0, col = "red")
p + stat_smooth()



################# the medium bugs version with locus, profile and dye effects and dose ##########################

bugsData = list(y = freq.df$obs, locus = freq.df$loc, dye = freq.df$dye, profile = freq.df$prof, N = length(freq.df$obs), X = freq.df$X,
                numLoci = 31, numProfiles = 102, numDyes = 4)
nChains = 1

bugsFile = here("medium.bugs.withDyeandDose.R")

## compile the model
system.time({sim = jags.model(file = bugsFile,
                              data = bugsData,
                              n.chains = nChains)})
update(sim, 10000)
parameters = c("Mu", "alpha.mu", "alpha.sigma", "alpha.locus", "gamma.dye", "mu", "beta.profile")
sim.sample = coda.samples(sim, parameters, n.iter = 1000)
simSummary = summary(sim.sample)

#graphing the locus effects
i = grep("^(alpha\\.locus).*$", rownames((simSummary$statistics)))
library(Hmisc)
errbar(1:31, simSummary$quantiles[i,3], simSummary$quantiles[i,1], simSummary$quantiles[i,5],
       xlab = "Locus", ylab = "Locus Effects")

#graphing the dye effects
i = grep("^(gamma).*$", rownames((simSummary$statistics)))
library(Hmisc)
numDyes = 4
errbar(1:numDyes, simSummary$quantiles[i,3], simSummary$quantiles[i,1], simSummary$quantiles[i,5],
       xlab = "Dye", ylab = "Dye Effects",
       axes = FALSE)
axis(2)
axis(1, at = 1:4)
box()

#graphing obs vs expected
i = grep("^(mu).*$", rownames((simSummary$statistics)))
fitted = simSummary$statistics[i,1] ## means
plot(fitted~bugsData$y, xlab = "Observed", ylab = "Fitted")
abline(c(0,1), col = "red")

results.df = data.frame(fitted = fitted, observed = bugsData$y)

library(ggplot2)
library(tidyverse)
p = results.df %>% ggplot(aes(x = observed, y = fitted)) + geom_point() + geom_abline(slope = 1, intercept = 0, col = "red")
p + stat_smooth()

fitted = simSummary$statistics[i,1] ## means
O = exp(bugsData$y)
E = exp(fitted)
plot(log10(O/E)~O, ylab = "log10(O/E)", xlab = "O")


