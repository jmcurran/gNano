#load required libraries
library(here)
library(rjags)

#loads the make bugs data method
source("makeBugsData.R")


################# formats the input data #########

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
                     obs = exp(obs[!is.na(obs)]),
                     dye = as.factor(dye[!is.na(obs)]),
                     X = X[!is.na(obs)])

## This model has a locus effect, and it treats the profile as a random effect
fit = aov(log(obs) ~ loc + dye + Error(prof) + X, data = freq.df)
summary(fit) ## and this shows we're well-justified to say there is a locus effect.
rm(obs)
rm(loc)
rm(prof)
rm(dye)
rm(X)



################# the medium bugs version with locus, profile and dye effects and dose ##########################

bugsData = list(y = freq.df$obs, N = length(freq.df$obs))

#, locus = freq.df$loc, dye = freq.df$dye, profile = freq.df$prof, N = length(freq.df$obs), X = freq.df$X,
#                numLoci = 31, numProfiles = 102, numDyes = 4)
nChains = 1

bugsFile = here("gamma.simple.R")

## compile the model
system.time({sim = jags.model(file = bugsFile,
                              data = bugsData,
                              n.chains = nChains)})
update(sim, 10000)
parameters = c("shape", "rate")
sim.sample = coda.samples(sim, parameters, n.iter = 1000)
simSummary = summary(sim.sample)


########## end of MCMC #### next section is all the graphing ##########
plot(density(bugsData$y, from = 0))
x = seq(0, 25000)
y = dgamma(x, shape = simSummary$statistics[2,1], rate = simSummary$statistics[1,1])
lines(x, y, col = "red")

################# locus effects ##########################

bugsData = list(y = freq.df$obs, N = length(freq.df$obs), locus = freq.df$loc,  numLoci = 31,
                pred = rep(NA, length(freq.df$obs)))
#dye = freq.df$dye, profile = freq.df$prof, N = length(freq.df$obs), X = freq.df$X,
#                numLoci = 31, numProfiles = 102, numDyes = 4)
nChains = 1

bugsFile = here("gamma.simple-locus.R")

## compile the model
system.time({sim = jags.model(file = bugsFile,
                              data = bugsData,
                              n.chains = nChains)})
system.time({update(sim, 10000)})
parameters = c("pred")
sim.sample = coda.samples(sim, parameters, n.iter = 1000)
simSummary = summary(sim.sample)

################# profile effects ##########################

bugsData = list(y = freq.df$obs, N = length(freq.df$obs), locus = freq.df$loc,  numLoci = 31,
                profile = freq.df$prof, numProfiles = 102,
                pred = rep(NA, length(freq.df$obs)))
#dye = freq.df$dye, profile = freq.df$prof, N = length(freq.df$obs), X = freq.df$X,
#                numLoci = 31, numProfiles = 102, numDyes = 4)
nChains = 1

bugsFile = here("gamma.simple-profile.R")

## compile the model
system.time({sim = jags.model(file = bugsFile,
                              data = bugsData,
                              n.chains = nChains)})
system.time({update(sim, 10000)})
parameters = c("pred", "Mu", "beta.profile")
sim.sample = coda.samples(sim, parameters, n.iter = 1000)
simSummary = summary(sim.sample)

b = grep("^beta.*$", rownames(simSummary$statistics))
beta.prof = simSummary$quantiles[b,]
library(Hmisc)
errbar(1:102, beta.prof[,3], beta.prof[,1], beta.prof[,5])
abline(col = "red", h = 0)

M = grep("Mu", rownames(simSummary$statistics))
simSummary$quantiles[M,]

pred = sim.sample[[1]][,-c(b,M)]

fit0 = glm(obs~1, data = freq.df, family = Gamma(link=log))
fit = glm(obs~as.factor(prof), data = freq.df, family = Gamma(link=log))


