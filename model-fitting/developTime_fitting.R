###############################################################################
# This code is part of the paper 
#
# Local thermal adaptation and local temperature regimes
# drive the performance of parasitic helminths under climate
# change: the case of Marshallagia marshalli from wild
# ungulates 
#
# by Aleuy et al. (2023, Global Change Biology)
#
# Author: Stephanie Peacock <stephanie.j.peacock at gmail.com>
# Data owner: O. Alejandro Aleuy <oaleuy@ucalgary.ca>
# March 4, 2022
#
###############################################################################


# File for loading data, exploring, and fitting models to estimate development and mortality
# of Marshallagia

# Load packages
library(here)
library(dclone)

# Load data
dat<-read.csv("data/Marsh_data.csv")

# Summary of data:
#		Location2: 6 locations
#		Temperature: up to 10 temperatures ( 5  7 10 15 20 25 30 33 35 38)
#		Incubator: 2 incubators for each Location2:Tempearture combination that has data
#		Well: Individual egg (one per well)
#		Hatching: The days since start of experiment that the egg hatched
#			Note: if Hatching==NA but L3_ocurrence!=NA, then Hatching_L3_BINO="YES" and eggs hatched as L3s (n=83)
# 	last_check_hatch: The days since the start of experiment that the egg was checked prior to hatching (t-dt)
# 	L3_ocurrence: The days since start of experiment that L3 occurred. 
#		last_check_L3: The days since the start of experiment that the larvae was checked prior to L3 (t-dt)
#		

# Create unique inclubator ID
incubatorLetter <- LETTERS[dat$Incubator]
dat$uniqueIncubator <- rep(NA, dim(dat)[1])
# Create exp number unique ID
dat$Expt <- rep(NA, dim(dat)[1])
	for(i in 1:dim(dat)[1]){
	if(dat$Location2[i] == "Colorado" | dat$Location2[i] == "Sheep_River_Colorado") dat$Expt[i] <- 1
	if(dat$Location2[i] == "Dawson" | dat$Location2[i] == "Sheep_River_Dawson") dat$Expt[i] <- 2
	if(dat$Location2[i] == "Sheep_Mountain" | dat$Location2[i] == "Sheep_River_ShM") dat$Expt[i] <- 3
	dat$uniqueIncubator[i] <- paste(dat$Temperature[i], "-", dat$Expt[i], incubatorLetter[i], sep="")
}
dat$numericIncubator <- as.numeric(factor(dat$uniqueIncubator))
# 56 unique inclubators

# Create true location data
dat$Location <- as.character(dat$Location2)
dat$Location[dat$Location == "Sheep_River_ShM" | dat$Location == "Sheep_River_Dawson" | dat$Location == "Sheep_River_Colorado"] <- "Sheep_River"

# Order factor from south to north
dat$Location <- factor(dat$Location, levels=c("Colorado", "Sheep_Mountain", "Sheep_River", "Dawson")) # Note this is wrong; Sheep_Mountain = South Yukon and Sheep_River = Alberta (!!)

# Sample size for each temperature/location combination
# n<-matrix(NA, nrow=length(unique(dat$Temperature)), ncol=length(unique(dat$Location2)))
# rownames(n) <- sort(unique(dat$Temperature))
# colnames(n) <- unique(dat$Location2)
# nj<-n
# 
# for(i in 1:dim(n)[1]){
# 	for (j in 1:dim(n)[2]){
# 		n[i,j]<-length(which(dat$Temperature==rownames(n)[i]&dat$Location2==unique(dat$Location2)[j]))
# 		nj[i,j]<-length(unique(dat$Incubator[which(dat$Temperature==rownames(n)[i]&dat$Location2==unique(dat$Location2)[j])]))
# }}

###############################################################################
# Development time to L3 
###############################################################################

# There are two options for how to treat individuals that survived and developed
# to L3 but had abnormal cells:
# # 1) treat these individuals as normal, surviving L3s (include.abnormal.as.dead == FALSE)
# # 2) treat these individuals as dead, with time of death as the time of 
#      development (include.abnormal.as.dead == TRUE)

include.abnormal.as.dead <- TRUE

dat$time.to.L3 <- dat$L3_ocurrence
dat$developmentCensored <- rep(0, dim(dat)[1])

#	If L3_ocurrence==NA, there could be several possibilities:
#		1) Dead_before_L3 - time that individual died, before they hatched to L3 -> treat as censored
	dat$developmentCensored[is.na(dat$Dead_before_L3) == FALSE] <- 1
	dat$time.to.L3[is.na(dat$Dead_before_L3) == FALSE] <- dat$Dead_before_L3[is.na(dat$Dead_before_L3) == FALSE]

#		2) Censored - individual survived to the end of the experiment but never hatched -> treat as censored
	dat$developmentCensored[which(dat$censored=="YES")] <- 1
	dat$time.to.L3[which(dat$censored=="YES")] <- dat$Final_day[which(dat$censored=="YES")]


dat.list<-list(
	nT = length(unique(dat$Temperature)),
	nL = 1,#length(unique(dat$Location)),
	nj = length(unique(dat$numericIncubator)),
	time = dat$time.to.L3,
	time_back = dat$last_check_L3,
	censored = dat$developmentCensored,
	temp = as.numeric(as.factor(dat$Temperature)),
	loc = rep(1, dim(dat)[1]),#as.numeric(dat$Location),
	incubator = dat$numericIncubator, #dat$Expt, 
	y = rep(1, length(dat$time.to.L3)),
	n = length(dat$time.to.L3)
	)

# To avoid errors, set NA time_back values to 1 
# Note: these do not come into model fitting as they are censored.
dat.list$time_back[dat.list$censored==1]<-1

# Add observed temperatures for MTE models
dat.list$T.obs <- sort(unique(dat$Temperature))

###############################################################################
# JAGS fitting
###############################################################################

# There are different models (independent temperature and MTE models)
# in the developTime_models.R file that are sourced below

source("developTime_models.R")

nChains <- 6 # January 23
nChains <- 8 # January 24
t.start <- numeric(5)
time2run <- numeric(5)

#------------------------------------------------------------------------------
# 1) Fit separate temperatures 
#------------------------------------------------------------------------------

cl <- makeCluster(nChains)
t.start[1]<-proc.time()[3]

fit <- jags.parfit(cl, data=dat.list, params=c("tau", "sigma_v", "sigma_e", "v"), n.adapt = 200000, n.update=100000, n.iter=5000, model=modelDevelop, n.chains=nChains)

stopCluster(cl)
time2run[1] <- round((proc.time()[3]-t.start[1])/(60), 1)
cat(paste("Process time =", time2run[1], "minutes"))

#------------------------------------------------------------------------------
# 2) Fit MTE - Arrhenius
#------------------------------------------------------------------------------
cl <- makeCluster(nChains)
t.start[2]<-proc.time()[3]

fitMTE.A <-jags.parfit(cl, data=dat.list, params=c("E", "sigma_e", "tau0", "sigma_v", "v"), model=modelDevelopMTE.A, n.adapt=200000, n.update=200000, n.iter=10000, n.chains=nChains) 

stopCluster(cl)
time2run[2] <- round((proc.time()[3]-t.start[2])/(60), 1)
cat(paste("Process time =", time2run[1], "minutes"))

# 82 mins with 200,000 - 100,000 - 5,000 but did not converge well

#------------------------------------------------------------------------------
# 3) Fit MTE - SSL - Jan 23 worked
#------------------------------------------------------------------------------
cl <- makeCluster(nChains)
t.start[3]<-proc.time()[3]

fitMTE.SSL<-jags.parfit(cl, data=dat.list, params=c("E", "El", "Tl", "sigma_e", "tau0", "sigma_v", "v"), model=modelDevelopMTE.SSL, n.adapt=200000, n.update=100000, n.iter=5000, n.chains=nChains)

stopCluster(cl)
time2run[3] <- round((proc.time()[3]-t.start[3])/(60), 1)
cat(paste("Process time =", time2run[5], "minutes"))
23 mins

#------------------------------------------------------------------------------
# 4) Fit MTE - SSU
# Doesn't make sense to fit this one because we know there's a lower bound.
# Just gives nonsensical values for activation energies etc.
#------------------------------------------------------------------------------
cl <- makeCluster(nChains)
t.start[4]<-proc.time()[3]

fitMTE.SSL<-jags.parfit(cl, data=dat.list, params=c("E", "Eh", "Th", "sigma_e", "tau0"), model=modelDevelopMTE.SSU, n.adapt=100000, n.update=50000, n.iter=1000, n.chains=nChains)

stopCluster(cl)
time2run[4] <- round((proc.time()[3]-t.start[4])/(60), 1)
cat(paste("Process time =", time2run[1], "minutes"))

#------------------------------------------------------------------------------
# 5) Fit MTE - SSUL
#------------------------------------------------------------------------------
cl <- makeCluster(nChains)
t.start[5]<-proc.time()[3]

fitMTE.SSUL <-jags.parfit(cl, data=dat.list, params=c("E", "El", "Eh", "Tl", "Th", "sigma_e", "tau0", "sigma_v", "v"), model = modelDevelopMTE.SSUL, n.adapt=200000, n.update=200000, n.iter = 10000, n.chains=nChains) 

stopCluster(cl)
time2run[5] <- round((proc.time()[3]-t.start[5])/(60), 1)
cat(paste("Process time =", time2run[5], "minutes"))

# 23 mins for 100,000 - 50,000 - 1000
# 147 mins for 200,000 - 100,000 - 50,000 with wrong tau0 (Jan 22)
# 69 mins for 100,000 - 50,000 - 5000 with corrected tau0 (works but not converged; Jan 23)
# 124 mins for 200,000 - 100,000 - 5000 but some convergence issues. More chains

#------------------------------------------------------------------------------
# Save output
#------------------------------------------------------------------------------

# save.image("ModelFits_SSUL_1loc_20190127.RData")
# Note: model output can be large and was not shared via GitHub.
# Email author for .RData files or iff there are issues running the above code.
