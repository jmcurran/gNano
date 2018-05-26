createGraphs <- function(sim.sample, saveDir, bugsData){

  #creates the directory if it doesn't already exist
  if (file.exists(saveDir)){
    #do nothing
  } else {
    #create it
    dir.create(file.path(saveDir))
  }
  
####### setting up the variables ##########
NumberLoci <- bugsData$numLoci
NumberDyes <- bugsData$numDyes
LocusNames <- c("K1", "M1", "R1", "K2", "Y1", "M2", "R2", "Y2", "Y3", "R3", "R4", "Y4", 
                "Y5", "R5", "Y6", "H1", "B1", "S1", "V1", "R6", "Y7", "Y8", "R7", "S2", 
                "R8", "Y9", "S3", "R9", "Y10", "R10", "Y11")
DyeNames <- c("A", "G", "C", "T")
DyeColours <- c("green", "blue", "black", "red")
alleles_at_locus <- bugsData$alleles_at_locus
profileData <- bugsData$P
profileDyes <- bugsData$profileDyes
NumberSamples <- dim(profileData)[[1]]
#sets the save dircetory (just so I can use old graphing code)
#saveDir = here("graphs//")



########pulls data out of jags model and into arrays for use in graphing code #########
lambda_capture <- NULL
try(
 lambda_capture <- sim.sample[[1]][,"lambda"]
)
#amp efficiency mean capture
A_mu_capture <- NULL
#counts the occurances of mus
#if(sum(colnames(sim.sample[[1]])=="mu.amp[1]")>0){
try(
  A_mu_capture <- sim.sample[[1]][,"mu.amp[1]"]
)
#}
#if(sum(colnames(sim.sample[[1]])=="mu.amp[2]")>0){
try(
 for (locus in 2:NumberLoci){
   A_mu_capture <- cbind(A_mu_capture, sim.sample[[1]][,paste("mu.amp[",locus,"]",sep="")])
 }
)
#}
#dye effect mean capture
D_mu_capture <- NULL
try(
 D_mu_capture <- sim.sample[[1]][,"mu.dye[1]"]
)
try(
 for (dye in 2:NumberDyes){
   D_mu_capture <- cbind(D_mu_capture, sim.sample[[1]][,paste("mu.dye[",dye,"]",sep="")])
 }
)
#amp efficiency sigma capture
A_sigma_capture <- NULL
try(
 A_sigma_capture <- sim.sample[[1]][,"sigma.sq.amp[1]"]
)
try(
 for (locus in 2:NumberLoci){
   A_sigma_capture <- cbind(A_sigma_capture, sim.sample[[1]][,paste("sigma.sq.amp[",locus,"]",sep="")])
 }
)
#dye effect sigma capture
D_sigma_capture <- NULL
try(
 D_sigma_capture <- sim.sample[[1]][,"sigma.sq.dye[1]"]
)
try(
 for (dye in 2:NumberDyes){
   D_sigma_capture <- cbind(D_sigma_capture, sim.sample[[1]][,paste("sigma.sq.dye[",dye,"]",sep="")])
 }
)
#expected peak hieghts
E <- array(NA,dim=c(NumberSamples, NumberLoci, 2))
for (sample in 1:NumberSamples){
  for (locus in 1:NumberLoci){
    if(!is.na(alleles_at_locus[sample,locus])){
      for (allele in 1:alleles_at_locus[sample,locus]){
        E[sample, locus, allele] <- mean(sim.sample[[1]][,paste("pred[", sample, ",", locus, ",", allele, "]",sep="")])
      }
    }
  }
}
T <- array(1,dim=c(NumberSamples))
for (sample in 1:NumberSamples){
  T[sample] <- mean(sim.sample[[1]][,paste("T[", sample, "]",sep="")])
}

#in case we want to cutr out the first part fo the graph
burnin <- 1
endcapture <- length(lambda_capture)





########## PLOT 3 - lambda over the MCMC ##########
if(!is.null(lambda_capture)){
ylim_max <- 1.1*max(lambda_capture)
plot(lambda_capture[1:endcapture], type="l", xlab="iteration (x50)", ylab="lambda", ylim=c(0,ylim_max))
dev.copy(jpeg, file=paste(saveDir, "lambda_MCMC.jpg", sep=""), height=2000, width=2000, res=300)
dev.off()


########## PLOT 4 - distribution of lambda after burnin ##########
binwidth <- hist(lambda_capture[burnin:endcapture])$mids[2]-hist(lambda_capture[burnin:endcapture])$mids[1]
hist(lambda_capture[burnin:endcapture], xlab="lambda", main=paste("lambda: N(", round(mean(lambda_capture[burnin:endcapture]),3), ",", round(sd(lambda_capture[burnin:endcapture]),3), ")",sep=""))
plot_x <- seq(min(lambda_capture[burnin:endcapture]), max(lambda_capture[burnin:endcapture]), length=1000)
lines(plot_x, binwidth*(endcapture - burnin)*dnorm(plot_x, mean(lambda_capture[burnin:endcapture]), sd(lambda_capture[burnin:endcapture])), col="red")
dev.copy(jpeg, file=paste(saveDir, "Lambda_histogram.jpg", sep=""), height=2000, width=2000, res=300)
dev.off()
}

########## PLOT 5 - amplification efficiencies means over the MCMC ##########
if(!is.null(A_mu_capture)){
ylim_max <- 1.1*max(A_mu_capture)
plot(A_mu_capture[1:endcapture,1], type="l", ylim=c(0, ylim_max), xlab="iteration", ylab="Locus amplification efficiency", xlim=c(1,(1.25*endcapture)))
for (locus in 2:NumberLoci){
  lines(A_mu_capture[1:endcapture,locus], ylim=c(0, ylim_max), xlim=c(1,(1.25*endcapture)))
}
#set up labels
x_label <- A_mu_capture[endcapture,]
y_label <- paste(LocusNames, " - ", signif(x_label,3),sep=" ")
legend("right", y_label[order(ordered(-x_label))], cex=0.7, y.intersp=0.5)

dev.copy(jpeg, file=paste(saveDir, "Amp_efficiencies_MCMC.jpg", sep=""), height=3000, width=3000, res=300)
dev.off()
}

########## PLOT 5.5 - amplification efficiencies sigmas over the MCMC ##########
if(!is.null(A_sigma_capture)){
plot(A_sigma_capture[1:endcapture,1], type="l", ylim=c(0.00000001, 4), xlab="iteration", ylab="Locus amplification efficiency variance", xlim=c(1,(1.25*endcapture)), log="y")
for (locus in 2:NumberLoci){
  lines(A_sigma_capture[1:endcapture,locus], ylim=c(0.00000001, 4), xlim=c(1,(1.25*endcapture)))
}
#set up labels
x_label <- A_sigma_capture[endcapture,]
y_label <- paste(LocusNames, " - ", signif(x_label,3),sep=" ")
legend("right", y_label[order(ordered(-x_label))], cex=0.7, y.intersp=0.5)
dev.copy(jpeg, file=paste(saveDir, "Amp_efficiencies_sigma_MCMC.jpg", sep=""), height=3000, width=4000, res=300)
dev.off()
}


########## PLOT 6 - dye effiencies means over the MCMC ##########
if(!is.null(D_mu_capture)){
plot(D_mu_capture[1:endcapture,1], type="l", ylim=c(0, 1.5), xlab="iteration", ylab="Dye amplification efficiency", xlim=c(1,(1.25*endcapture)))
for (locus in 2:NumberDyes){
  lines(D_mu_capture[1:endcapture,locus], ylim=c(0, 1.5), xlim=c(1,(1.25*endcapture)))
}
#set up labels
x_label <- D_mu_capture[endcapture,]
y_label <- paste(DyeNames, " - ", signif(x_label,3),sep=" ")
legend("right", y_label[order(ordered(-x_label))])
dev.copy(jpeg, file=paste(saveDir, "Dye_effects_MCMC.jpg", sep=""), height=2000, width=3000, res=300)
dev.off()
}

########## PLOT 6.5 - dye effiencies variances over the MCMC ##########
if(!is.null(D_sigma_capture)){
plot(D_sigma_capture[1:endcapture,1], type="l", ylim=c(0, 0.05), xlab="iteration", ylab="Dye amplification efficiency variance", xlim=c(1,(1.25*endcapture)))
for (locus in 2:NumberDyes){
  lines(D_sigma_capture[1:endcapture,locus], ylim=c(0, 0.05), xlim=c(1,(1.25*endcapture)))
}
#set up labels
x_label <- D_sigma_capture[endcapture,]
y_label <- paste(DyeNames, " - ", signif(x_label,3),sep=" ")
legend("right", y_label[order(ordered(-x_label))])
dev.copy(jpeg, file=paste(saveDir, "Dye_effects_sigma_MCMC.jpg", sep=""), height=2000, width=3000, res=300)
dev.off()
}

########## PLOT 7 - locus amplification efficiency distributions##########
if(!is.null(A_mu_capture) && !is.null(A_sigma_capture)){
simulated_curves <- 1000
plot_x <- seq(0.01, 10, length=1000)
color=rgb(0,0,0,alpha=0.1)
for (locus in 1:NumberLoci){
  #calculates mean and varaince of aplification efficiency mean
  A_mu_mean <- mean(A_mu_capture[burnin:endcapture,locus])
  A_mu_sigma <- sd(A_mu_capture[burnin:endcapture,locus])
  #calculates mean and variance of aplification efficiency variance
  A_sigma_mean <- mean(A_sigma_capture[burnin:endcapture,locus])
  A_sigma_sigma <- sd(A_sigma_capture[burnin:endcapture,locus])
  ylim_max <- 3*dlnorm(A_mu_mean, log(A_mu_mean), sqrt(A_sigma_mean))
  #opens blank plot
  plot(0,0, xlim=c(0, (2*A_mu_mean)), ylim=c(0,ylim_max), xlab="Amplification effiency", ylab="", main = paste("locus ", LocusNames[locus], " amplification efficiency", sep=""))
  #randomly samples from mean and variance distributions 
  for (randcurve in 1:simulated_curves){
    random_A_mu <- -1
    while (random_A_mu <= 0){
      random_A_mu <- rnorm(1, A_mu_mean, sqrt(A_mu_sigma))
    }
    random_A_sigma <- -1
    while(random_A_sigma <= 0){
      random_A_sigma <- rnorm(1, A_sigma_mean, sqrt(A_sigma_sigma))
    }
    #draw LN
    lines(plot_x, dlnorm(plot_x, log(random_A_mu), sqrt(random_A_sigma)), col=color, xlim=c(0, (2*A_mu_mean)), ylim=c(0,ylim_max))
  }
  #plots mean values
  lines(plot_x, dlnorm(plot_x, log(A_mu_mean), sqrt(A_sigma_mean)), col="red", xlim=c(0, (2*A_mu_mean)), ylim=c(0,ylim_max))
  #saves plots
  dev.copy(jpeg, file=paste(saveDir, LocusNames[locus], " Amp efficiency.jpg", sep=""), height=2000, width=2000, res=300)
  dev.off()
}
}

########## PLOT 8 - dye amplification efficiency distributions##########
if(!is.null(D_mu_capture) && !is.null(D_sigma_capture)){
simulated_curves <- 1000
plot_x <- seq(0.01, 10, length=1000)
color=rgb(0,0,0,alpha=0.1)
#opens blank plot
ylim_max <- 10
plot(-1000,-1000, xlim=c(0, 2), ylim=c(0,ylim_max), xlab="Dye amplification effiency", ylab="", main = "")
for (dye in 1:NumberDyes){
  #color <- col2rgb(DyeColours[dye], alpha=0.01)
  #calculates mean and varaince of aplification efficiency mean
  D_mu_mean <- mean(D_mu_capture[burnin:endcapture,dye])
  D_mu_sigma <- sd(D_mu_capture[burnin:endcapture,dye])
  #calculates mean and variance of aplification efficiency variance
  D_sigma_mean <- mean(D_sigma_capture[burnin:endcapture,dye])
  D_sigma_sigma <- sd(D_sigma_capture[burnin:endcapture,dye])
  #randomly samples from mean and variance distributions 
  for (randcurve in 1:simulated_curves){
    random_D_mu <- -1
    while(random_D_mu <= 0){
      random_D_mu <- rnorm(1, D_mu_mean, sqrt(D_mu_sigma))
    }
    random_D_sigma <- -1
    while(random_D_sigma <= 0){
      random_D_sigma <- rnorm(1, D_sigma_mean, sqrt(D_sigma_sigma))
    }
    #draw LN
    lines(plot_x, dlnorm(plot_x, log(random_D_mu), sqrt(random_D_sigma)), col=color, xlim=c(0, 2), ylim=c(0,ylim_max))
  }
  #draw curve using mean value
  lines(plot_x, dlnorm(plot_x, log(D_mu_mean), sqrt(D_sigma_mean)), col=DyeColours[dye], xlim=c(0, 2), ylim=c(0,ylim_max))
}
#saves plots
dev.copy(jpeg, file=paste(saveDir, "Dye efficiency.jpg", sep=""), height=2000, width=2000, res=300)
dev.off()
}


########## PLOT 9 - Expected peak heights ##########
#modelling peak heights
mean_peak_height <- mean(profileData, na.rm=TRUE)
peak_count <- sum(profileData>0, na.rm=TRUE)
binwidth <- hist(profileData)$mids[2]-hist(profileData)$mids[1]
hist(profileData, main="observed peak heights", xlab="peak height (rfu)", xlim=c(0,30000))
plot_x <- seq(1,30000, length = 10000)
plot_y <- dexp(plot_x, 1/mean_peak_height)
lines(plot_x, binwidth*peak_count*plot_y, col="red")
dev.copy(jpeg, file=paste(saveDir, "Observed_peak_heights.jpg", sep=""), height=1500, width=3000, res=300)
dev.off()


###################   Model checks   #################



########## PLOT 10 - expected and observed peak height variability ##########
E_over_O <- E/profileData
xlim_range <- round(max(T), -3) #rounded to thousands
plot(-1000,-1000, xlim=c(0, xlim_range), ylim=c(0,4), xlab="Template", ylab="O/E", main = paste("Peak Height varibility", sep=""))
for (sample in 1:NumberSamples){
  points(seq(T[sample], length=length(E_over_O[sample,,])), E_over_O[sample,,], xlim=c(0, xlim_range), ylim=c(0,4))
}
#now draw lines from analysis
mean_lambda <- mean(lambda_capture[burnin:endcapture])
plot_seq <- seq(1, xlim_range, length=xlim_range)
#percentage to capture
capture_percentage <- 95
Z_value <- qnorm((1 - capture_percentage/100)/2, lower.tail=FALSE)
lines(plot_seq , exp(sqrt(Z_value*mean_lambda/plot_seq)), col="red")
lines(plot_seq , exp(-sqrt(Z_value*mean_lambda/plot_seq)), col="red")
#determine actual coverage
count_within_bounds <- 0
count_outside_bounds <- 0
for (sample in 1:NumberSamples){
  for (locus in 1:NumberLoci){
    if(!is.na(alleles_at_locus[sample,locus])){
    for (allele in 1:alleles_at_locus[sample,locus]){
      if(!is.na(E_over_O[sample, locus, allele])) {
        if (E_over_O[sample, locus, allele] <= exp(sqrt(Z_value*mean_lambda/T[sample])) && E_over_O[sample, locus, allele] >= exp(-sqrt(Z_value*mean_lambda/T[sample]))){
          count_within_bounds <- count_within_bounds +1
        }
        else {
          count_outside_bounds <- count_outside_bounds + 1
        }
      }
    }
    }
  }
}
percent_within <- round(100*count_within_bounds/(count_within_bounds+count_outside_bounds), 2)
legend("topright", c("O/E", paste(capture_percentage, " % bounds (",percent_within,"% actual)",sep="")), pch=c("o","-"), col=c("black", "red"))
dev.copy(jpeg, file=paste(saveDir, "Peak_Height_Variability.jpg", sep=""), height=2000, width=3000, res=300)
dev.off()


########## PLOT 11 - expected and observed Heterozygote balance by dye ##########
if(!is.null(D_mu_capture) && !is.null(D_sigma_capture)){
#divide first peak by second in all profiles and for all loci
Hb <- profileData[,,1]/profileData[,,2]
#need 6 box plots  A/C, A/G, A/T, C/G, C/T, G/T
#number equivalent 1/3, 1/2, 1/4, 3/2, 3/4, 2/4
combinations <- 6
dye1 <- c(1, 1, 1, 3, 3, 2)
dye2 <- c(3, 2, 4, 2, 4, 4)
#now draws simulated vs expected ratios
for (combo in 1:combinations){
  dye1_dye2_observed <- Hb[(profileDyes[,,1]==dye1[combo] & profileDyes[,,2]==dye2[combo])]
  dye1_mean <- mean(D_mu_capture[burnin:endcapture,dye1[combo]])
  dye1_sigma <- sd(D_sigma_capture[burnin:endcapture,dye1[combo]])
  dye1_simulated <- rlnorm(1000, log(dye1_mean), sqrt(dye1_sigma))
  dye2_mean <- mean(D_mu_capture[burnin:endcapture,dye2[combo]])
  dye2_sigma <- sd(D_sigma_capture[burnin:endcapture,dye2[combo]])
  dye2_simulated <- rlnorm(1000, log(dye2_mean), sqrt(dye2_sigma))
  #simulates peak height variability from stochastic effects (random peak between 500 and 5000 rfu)
  peak1_variability <- array(1,dim=c(1000))
  for (var in 1:1000){
    #random peak height from observed peak heights, which are exponential, but above AT
    peak <- max(50, rexp(1, 1/mean_peak_height))
    peak1_variability[var] <- rlnorm(1, log(peak), sqrt(mean_lambda/peak))/rlnorm(1, log(peak), sqrt(mean_lambda/peak))
  }
  dye1_dye2_simulated <- peak1_variability*dye1_simulated/dye2_simulated
  if (combo == 1){
    boxplot(dye1_dye2_observed, at=1, xlim=c(1,combinations*2), ylim=c(0,7), main="observed (O) and simulated (S) dye balance", xlab="", ylab="Ratio of Heights", col="red")
  }
  else { #just add
    boxplot(dye1_dye2_observed, at=((combo-1)*2+1), xlim=c(1,combinations*2), ylim=c(0,7), col="red", add=TRUE)
  }
  boxplot(dye1_dye2_simulated, at=(combo*2), xlim=c(1,combinations*2), ylim=c(0,7), col="blue", add=TRUE)
}
axis(1, at=1:(2*combinations), labels=c("A/C(O)", "A/C(S)", "A/G(O)", "A/G(S)", "A/T(O)", "A/T(S)", "C/G(O)", "C/G(S)", "C/T(O)", "C/T(S)", "G/T(O)", "G/T(S)"))
dev.copy(jpeg, file=paste(saveDir, "Dye_amplification_Efficiency.jpg", sep=""), height=2000, width=6000, res=300)
dev.off()
}


########## PLOT 12 - Peak height per locus (relatve to avergae in profile) ##########
if(!is.null(A_mu_capture) && !is.null(A_sigma_capture)){
Av_peakheight_per_profile <- apply(profileData, c(1), mean, na.rm=TRUE)
scaled_O  <- array(data=NA, dim=c(NumberSamples, NumberLoci, 2))
#add across peaks of heterozygote
scaled_O  <- apply(profileData, c(1, 2), sum, na.rm=TRUE)
for (profile in 1:NumberSamples){
  scaled_O[profile,] <- scaled_O[profile,]/Av_peakheight_per_profile[profile]
}
#plot as boxplot
for (locus in 1:NumberLoci){
  if (locus == 1){
    boxplot(scaled_O[,locus], at=(2*(locus-1)), xlim=c(1,NumberLoci*2), ylim=c(0,6), main="observed (O) and simulated (S) lous amp balance", xlab="", ylab="ratio peak heights / average", col="red")
  }
  else {
    boxplot(scaled_O[,locus], at=(2*(locus-1)), xlim=c(1,NumberLoci*2), ylim=c(0,6), col="red", add=TRUE)
  }
}
#now the simulated values
random_profile <- array(0, dim=c(1000, NumberLoci))
for (var in 1:1000){
  #choose template
  random_template <- max(50, rexp(1, 1/mean_peak_height))
  for (locus in 1:NumberLoci){
    A_mu_means <- mean(A_mu_capture[burnin:endcapture,locus])
    A_sigma_means <- mean(A_sigma_capture[burnin:endcapture,locus])
    #random peak height variability value
    peak_variability <- rlnorm(1, log(random_template), sqrt(mean_lambda/random_template))
    #random locus amp efficiency value
    random_locus_amp <- rlnorm(1, log(A_mu_means), sqrt(A_sigma_means))
    random_profile[var,locus] <- peak_variability*random_locus_amp
  }
  #scale
  random_profile[var,] <- random_profile[var,]/mean(random_profile[var,])
}
#plot simulated data
for (locus in 1:NumberLoci){
  boxplot(random_profile[,locus], at=(2*(locus-1)+1), xlim=c(1,NumberLoci*2), ylim=c(0,6), col="blue", add=TRUE)
}
axis(1, at=seq(1,(2*NumberLoci-1), length = NumberLoci), labels=LocusNames)
dev.copy(jpeg, file=paste(saveDir, "Locus_amplification_Efficiency.jpg", sep=""), height=2000, width=6000, res=300)
dev.off()
}


#GR converence on diagnostics
library(coda)


########## PLOT 14 - Gelman Rubin convergence diagnostics for other model parameters ##########
#values for graph
number_points <- 1 + 1 + NumberDyes + NumberDyes + NumberLoci + NumberLoci
graph_values <- array(data=0, dim=c(number_points))
graph_labels <- array(data=0, dim=c(number_points))
graph_counter <- 1
#dataset probability (not applicable)
graph_values[graph_counter] <- 0#GR_p_total
graph_labels[graph_counter] <- "dataset probability"
graph_counter <- graph_counter + 1
#lambda
GR <- gelman.diag(sim.sample[,"lambda"])
GR_lambda <- GR[[1]][1]
graph_values[graph_counter] <- GR_lambda
graph_labels[graph_counter] <- "lambda"
graph_counter <- graph_counter + 1
#D mu
GR_D_mu <- array(data=NA, dim=c(NumberDyes))
for (dye in 1:NumberDyes){
  try(
    GR <- gelman.diag(sim.sample[,paste("mu.dye[",dye,"]",sep="")])
  )
  try(
    GR_D_mu[dye] <- GR[[1]][1]
  )
  try(
    graph_values[graph_counter] <- GR_D_mu[dye]
  )
  graph_labels[graph_counter] <- paste("Dye amp mean (", DyeNames[dye], ")", sep="")
  graph_counter <- graph_counter + 1
}
#D sigma
GR_D_sigma <- array(data=NA, dim=c(NumberDyes))
for (dye in 1:NumberDyes){
  try(
   GR <- gelman.diag(sim.sample[,paste("sigma.sq.dye[",dye,"]",sep="")])
  )
  try(
   GR_D_sigma[dye] <- GR[[1]][1]
  )
  try(
   graph_values[graph_counter] <- GR_D_sigma[dye]
  )
  graph_labels[graph_counter] <- paste("Dye amp var (", DyeNames[dye], ")", sep="")
  graph_counter <- graph_counter + 1
}
#A mu
GR_A_mu <- array(data=NA, dim=c(NumberLoci))
for (locus in 1:NumberLoci){
  try(
  GR <- gelman.diag(sim.sample[,paste("mu.amp[",dye,"]",sep="")])
  )
  try(
  GR_A_mu[locus] <- GR[[1]][1]
  )
  try(
  graph_values[graph_counter] <- GR_A_mu[locus]
  )
  graph_labels[graph_counter] <- paste("Locus amp mean (", LocusNames[locus], ")", sep="")
  graph_counter <- graph_counter + 1
}
#A sigma
GR_A_sigma <- array(data=NA, dim=c(NumberLoci))
for (locus in 1:NumberLoci){
  try(
  GR <- gelman.diag(sim.sample[,paste("sigma.sq.amp[",dye,"]",sep="")])
  )
  try(
  GR_A_sigma[locus] <- GR[[1]][1]
  )
  try(
  graph_values[graph_counter] <- GR_A_sigma[locus]
  )
  graph_labels[graph_counter] <- paste("Locus amp var (", LocusNames[locus], ")", sep="")
  graph_counter <- graph_counter + 1
}
graph_labels <- paste(graph_labels, " : ", round(graph_values,2), sep="")
#now draw the graph
plot (graph_values, main="Gelman Rubin convergence diagnostic", xlab="Parameter", ylab="GR value", xaxt='n', ylim=c(0,12))
abline(h=1.2, col="grey")
text(seq(1,(graph_counter-1), length=(graph_counter-1)), (graph_values+2.5), labels=graph_labels, pos=3, cex=0.7, srt=90)
dev.copy(jpeg, file=paste(saveDir, "parameter_GR.jpg", sep=""), height=2000, width=6000, res=300)
dev.off()

#################### END PLOTS ################################


} #end function
