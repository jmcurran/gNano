source("run_Gnano_model.R")

#the uniform full model
run_Gnano_model("gnano_truncated.bugs.R", "uniform")

#the uniform full model modified dye prior
run_Gnano_model("gnano_truncated_uniform_mod.bugs.R", "uniform_mod")


#the gamma full model
run_Gnano_model("gnano_truncated-gamma.bugs.R", "fullGamma")

#the gamma noAmp model
run_Gnano_model("gnano_truncated-gamma.no.amp.bugs.R", "fullGammaNoAmp")

#the gamma noDye model
run_Gnano_model("gnano_truncated-gamma.no.dye.bugs.R", "fullGammaNoDye")

#the gamma noAmpNoDye model
run_Gnano_model("gnano_truncated-gamma.no.amp.no.dye.bugs.R", "fullGammaNoAmpNoDye")

#the gamma full gammaPredictor model
run_Gnano_model("gnano_truncated-gamma.gammaPred.bugs.R", "fullGammaPredGamma")

#the gamma full constantAmpSigma model
run_Gnano_model("gnano_truncated-gamma.OneAmpSig.bugs.R", "fullGammaConstantAmpSigma")

#the gamma full constantDyeSigma model
run_Gnano_model("gnano_truncated-gamma.OneDyeSig.bugs.R", "fullGammaConstantDyeSigma")

#the gamma full constantAmpAndDyeSigma model
run_Gnano_model("gnano_truncated-gamma.OneAmpAndDyeSig.bugs.R", "fullGammaConstantAmpAndDyeSigma")

