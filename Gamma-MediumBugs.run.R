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

bugsData = list(y = freq.df$obs, locus = freq.df$loc, dye = freq.df$dye, profile = freq.df$prof, N = length(freq.df$obs), X = freq.df$X,
                numLoci = 31, numProfiles = 102, numDyes = 4)
nChains = 1

bugsFile = here("gamma.medium.bugs.withDyeandDose.R")

## compile the model
system.time({sim = jags.model(file = bugsFile,
                              data = bugsData,
                              n.chains = nChains)})
update(sim, 10000)
parameters = c("Mu", "alpha.mu", "alpha.sigma", "alpha.locus", "gamma.dye", "mu", "beta.profile", "tau")
sim.sample = coda.samples(sim, parameters, n.iter = 1000)
simSummary = summary(sim.sample)


########## end of MCMC #### next section is all the graphing ##########


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
plot(log(fitted) - log(bugsData$y) ~ log(fitted), xlab = expression(hat(y)), ylab = expression(hat(y)-y))
h = abline(h = 0, col = "red")

results.df = data.frame(fitted = fitted, observed = bugsData$y)
results.df = results.df %>% 
  mutate(residuals = fitted - observed)
p = results.df %>% ggplot(aes(x = fitted, y = abs(residuals))) + geom_point() + stat_smooth()


library(ggplot2)
library(tidyverse)
p = results.df %>% ggplot(aes(x = observed, y = fitted)) + geom_point() + geom_abline(slope = 1, intercept = 0, col = "red")
p + stat_smooth()

fitted = simSummary$statistics[i,1] ## means
#getting the template values
T <- NULL
E <- NULL
for (i in 1:length(freq.df$obs)){
  T <- c(T, mean(sim.sample[[1]][, paste("beta.profile[", as.numeric(freq.df$prof[i]), "]", sep ="")]))
  E <- c(E, mean(sim.sample[[1]][, paste("mu[", i, "]", sep ="")]))
}
T <- exp(T)
#observed
O = exp(bugsData$y)
#expected
E = exp(E)
#plot(log10(O/E)~E, ylab = "log10(O/E)", xlab = "E")
plot(log10(O/E)~T, ylab = "log10(O/E)", xlab = "T")
abline(h=0, col="red")
