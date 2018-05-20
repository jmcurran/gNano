## run after bugsData has been set up

## the full model with uniform priors
bugsFile = here("gnano_truncated.bugs.R")
#compile the model
sim = jags.model(file = bugsFile, data = bugsData, n.chain = 4)
dic.unif <- dic.samples(model = sim, n.iter = 50000, thin = 50, type = "pD")

## the full model with inverse gamma priors
bugsFileGamma = here("gnano_truncated-gamma.bugs.R")
#compile the model
simGamma = jags.model(file = bugsFileGamma, data = bugsData, n.chain = 4)
dic.gamma <- dic.samples(model = simGamma, n.iter = 50000, thin = 50, type = "pD")

## the no amp inverse gamma priors model
bugsFileGammaNoAmp = here("gnano_truncated-gamma.no.amp.bugs.R")
#compile the model
simGammaNoAmp = jags.model(file = bugsFileGammaNoAmp, data = bugsData, n.chain = 4)
dic.gamma.NoAmp <- dic.samples(model = simGammaNoAmp, n.iter = 50000, thin = 50, type = "pD")

## the no dye inverse gamma priors model
bugsFileGammaNoDye = here("gnano_truncated-gamma.no.dye.bugs.R")
#compile the model
simGammaNoDye = jags.model(file = bugsFileGammaNoDye, data = bugsData, n.chain = 4)
dic.gamma.NoDye <- dic.samples(model = simGammaNoDye, n.iter = 50000, thin = 50, type = "pD")

## the no amp and no dye inverse gamma priors model
bugsFileGammaNoAmpNoDye = here("gnano_truncated-gamma.no.amp.no.dye.bugs.R")
#compile the model
simGammaNoAmpNoDye = jags.model(file = bugsFileGammaNoAmpNoDye, data = bugsData, n.chain = 4)
dic.gamma.NoAmpNoDye <- dic.samples(model = simGammaNoAmpNoDye, n.iter = 50000, thin = 50, type = "pD")

## printing out the DIC results
print(c("standard uniform model:"))
dic.unif
print(c("standard gamma model:"))
dic.gamma
print(c("No amp gamma model:"))
dic.gamma.NoAmp
print(c("No dye gamma model:"))
dic.gamma.NoDye
print(c("No amp and no dye gamma model:"))
dic.gamma.NoAmpNoDye

