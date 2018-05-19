## run after both uniform and gamma bugs setups

## model comparison code
simGamma = jags.model(file = bugsFileGamma, data = bugsData, n.chain = 4)
sim = jags.model(file = bugsFile, data = bugsData, n.chain = 4)

## now do model comparison using DIC
dic.unif <- dic.samples(model = sim, n.iter = 50000, thin = 50, type = "pD")
dic.gamma <- dic.samples(model = simGamma, n.iter = 50000, thin = 50, type = "pD")
DIC.comparison <- diffdic(dic.unif, dic.gamma)
DIC.comparison


