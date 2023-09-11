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
dat$Location <- factor(dat$Location, levels=c("Colorado", "Sheep_River", "Sheep_Mountain", "Dawson"))

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

hist(dat$Dead_before_L3)
hist(dat$last_check_Mortality)

# eggs that hatched
sum(is.na(dat$Hatching)==FALSE)
	# and developed to L3
	sum(is.na(dat$Hatching)==FALSE & is.na(dat$L3_ocurrence) == FALSE)
	# died before L3
	sum(is.na(dat$Hatching)==FALSE & is.na(dat$L3_ocurrence) == TRUE & is.na(dat$Dead_before_L3) == FALSE)
	# survived to end of expt but never developed to L3
	sum(is.na(dat$Hatching)==FALSE & is.na(dat$L3_ocurrence) == TRUE & is.na(dat$Dead_before_L3) == TRUE)
	# developed to L3 but abnormal
	sum(is.na(dat$L3_ocurrence) == FALSE & dat$abnormal_cells == "YES", na.rm=TRUE)

# eggs that hatched as L3s
length(which(is.na(dat$Hatching)==TRUE & is.na(dat$L3_ocurrence) == FALSE))
#***** row 349: Colorado 5C, Incubator 2, Well d3 has Hatching = NA and L3_ocurrence = 249 , but Hatching_L3_BINO = "NO"
length(which(dat$Hatching_L3_BINO == "YES"))

# Eggs that never hatched
dat$notHatched <- rep(0, dim(dat)[1])
dat$notHatched[which(is.na(dat$Hatching)==TRUE & is.na(dat$L3_ocurrence) == TRUE)] <- 1
sum(dat$notHatched)
	# Died as eggs
	length(which(dat$notHatched == 1 & is.na(dat$Dead_before_L3) == FALSE))
	head(dat[which(dat$notHatched == 1 & is.na(dat$Dead_before_L3) == FALSE),])
	# Survived to end as eggs
	length(which(dat$notHatched == 1 & is.na(dat$Dead_before_L3) == TRUE))
	
head(dat[which(is.na(dat$L3_ocurrence)==TRUE & is.na(dat$Dead_before_L3)==TRUE & is.na(dat$Hatching)==TRUE),])


###############################################################################
# Survival to L3 
###############################################################################

# There are two options for how to treat individuals that survived and developed
# to L3 but had abnormal cells:
# # 1) treat these individuals as normal, surviving L3s (include.abnormal.as.dead == FALSE)
# # 2) treat these individuals as dead, with time of death as the time of 
#      development (include.abnormal.as.dead == TRUE)

include.abnormal.as.dead <- TRUE


if(include.abnormal.as.dead == FALSE) {
	#	Option #1If Dead_before_L3 == NA, the data are considered censored.
	dat$deathCensored <- rep(0, dim(dat)[1])
	dat$deathCensored[is.na(dat$Dead_before_L3)] <- 1
	cat(paste("Number of censored individuals:", sum(dat$deathCensored))) # 1667 censored
	
	# Censoring can be do to either
	#		1) the individual survived to L3 and therefore L3_occurrence != NA. In this case, 
	#			the censored time is the time that they developed.
	C1 <- which(dat$deathCensored == 1 & is.na(dat$L3_ocurrence) == FALSE)
	dat$last_check_Mortality[C1] <- dat$L3_ocurrence[C1]
	
	#		2) the inidividual survived to the end of the experiment but never developed to L3
	# 		 (Of those 98 individuals in category 2, there were 36 who never even hatched.)
	C2 <- which(dat$deathCensored == 1 & is.na(dat$L3_ocurrence) == TRUE)
	dat$last_check_Mortality[C2] <- dat$Final_day[C2]
	
} else if (include.abnormal.as.dead == TRUE) {
	# Option #2: If Dead_before_L3 == NA & not abnormal, then censored
	# BUT if L3 cells were abnormal, then consider time of death to be
	# the L3_occurrence
	dat$deathCensored <- rep(0, dim(dat)[1])
	dat$deathCensored[is.na(dat$Dead_before_L3) & (dat$abnormal_cells == "NO" | is.na(dat$abnormal_cells) == TRUE)] <- 1
	cat(paste("Number of censored individuals:", sum(dat$deathCensored))) # 1555 censored
	
	# Censoring can be do to either
	#		1) the individual survived to L3 AND was normal, and therefore L3_occurrence != NA. In this case, 
	#			the censored time is the time that they developed.
	C1 <- which(dat$deathCensored == 1 & is.na(dat$L3_ocurrence) == FALSE & dat$abnormal_cells == "NO")
	dat$last_check_Mortality[C1] <- dat$L3_ocurrence[C1]
	
	#		2) the inidividual survived to the end of the experiment but never developed to L3
	# 		 (Of those 98 individuals in category 2, there were 36 who never even hatched.)
	C2 <- which(dat$deathCensored == 1 & is.na(dat$L3_ocurrence) == TRUE)
	dat$last_check_Mortality[C2] <- dat$Final_day[C2]
	
	# 	3) If individual survived to L3 but was abnormal, set time of death as time of development
	#			 to L3 (uncensored)
	UC1 <- which(dat$abnormal_cells == "YES")
	dat$last_check_Mortality[UC1] <- dat$last_check_L3[UC1]
	dat$Dead_before_L3[UC1] <- dat$L3_ocurrence[UC1]
}



#------------------------------------------------------------------------------
# Prepare data list for JAGS
#------------------------------------------------------------------------------

dat.list<-list(
	nT = length(unique(dat$Temperature)),
	nL = length(unique(dat$Location)),
	nj = length(unique(dat$numericIncubator)),
	time = dat$Dead_before_L3,
	time_back = dat$last_check_Mortality,
	censored = dat$deathCensored,
	temp = as.numeric(as.factor(dat$Temperature)),
	loc = as.numeric(dat$Location),#rep(1, dim(dat)[1])
	incubator = dat$numericIncubator, #dat$Expt, 
	y = rep(1, length(dat$Dead_before_L3)),
	n = length(dat$Dead_before_L3)
	)

# Check that all censored individuals have time as NA
sum(is.na(dat.list$time[dat.list$censored == 1])) == sum(dat.list$censored)
# but all individuals have a time_back
sum(is.na(dat.list$time_back))

# Replace time for censored with 1 to avoid errors (doesn't come into likelihood calc)
dat.list$time[dat.list$censored == 1] <- 1

# Add observed temperatures for MTE models
dat.list$T.obs <- sort(unique(dat$Temperature))

meanTime2Death <- matrix(NA, nrow=dat.list$nT, ncol=dat.list$nL)
nTime2Death <- matrix(NA, nrow=dat.list$nT, ncol=dat.list$nL)
for(i in 1:dat.list$nT){
	for(j in 1:dat.list$nL){
		meanTime2Death[i,j] <- mean(dat.list$time_back[dat.list$temp == i & dat.list$loc == j], na.rm=TRUE)
		nTime2Death[i,j] <- sum(dat.list$censored[dat.list$temp == i & dat.list$loc == j] == 0)
	}}

cols <- c(Col = '#e41a1c', ShR = '#377eb8', ShM = '#4daf4a', Daw = '#984ea3')
bp <- barplot(t(meanTime2Death), beside=TRUE, ylab="Time to death (days)", xlab=expression(paste("Temperature (", degree, "C)")), col=cols, ci.width = 0.2, las=1)
text(bp, t(meanTime2Death), t(nTime2Death), cex=0.6, pos=3, font=2, xpd=NA)
abline(h = 0)
legend('topright', fill=cols, legend=levels(dat$Location), ncol=1, bty="n")

# #------------------------------------------------------------------------------
# # Approach #2: survival through life stage (D_L in Molnar et al. 2013)
# 
# dat.list<-list(
# 	nT = length(unique(dat$Temperature)),
# 	nL = length(unique(dat$Location)),
# 	nj = length(unique(dat$numericIncubator)),
# 	time = dat$Dead_before_L3,
# 	time_back = dat$last_check_Mortality,
# 	censored = dat$deathCensored,
# 	temp = as.numeric(as.factor(dat$Temperature)),
# 	loc = as.numeric(dat$Location),#rep(1, dim(dat)[1])
# 	incubator = dat$numericIncubator, #dat$Expt, 
# 	y = rep(1, length(dat$Dead_before_L3)),
# 	n = length(dat$Dead_before_L3)
# )
#------------------------------------------------------------------------------
# JAGS fitting
#------------------------------------------------------------------------------

# Source file with all five models that were fit
source("surv_models.R")

nChains <- 8 # January 24
t.start <- numeric(5)
time2run <- numeric(5)

#------------------------------------------------------------------------------
# 1) Fit separate temperatures
#------------------------------------------------------------------------------

cl <- makeCluster(nChains)
t.start[1]<-proc.time()[3]

fit <- jags.parfit(cl, data=dat.list, params=c("mu"), n.adapt=5000, n.update=5000, n.iter=5000, model=modelSurv, n.chains=nChains)

stopCluster(cl)
time2run[1] <- round((proc.time()[3]-t.start[1])/(60), 1)
cat(paste("Process time =", time2run[1], "minutes"))

#------------------------------------------------------------------------------
# 2) Fit MTE - Arrhenius
#------------------------------------------------------------------------------
cl <- makeCluster(nChains)
t.start[2]<-proc.time()[3]

fitMTE.A<-jags.parfit(cl, data=dat.list, params=c("E", "mu0"), model=modelSurvMTE.A, n.adapt=5000, n.update=5000, n.iter=5000, n.chains=nChains) 

stopCluster(cl)
time2run[2] <- round((proc.time()[3]-t.start[2])/(60), 1)
cat(paste("Process time =", time2run[1], "minutes"))

#------------------------------------------------------------------------------
# 3) Fit MTE - SSU
#------------------------------------------------------------------------------
cl <- makeCluster(nChains)
t.start[3]<-proc.time()[3]

fitMTE.SSU <-jags.parfit(cl, data=dat.list, params=c("E", "Eh", "Th", "mu0"), model=modelSurvMTE.SSU, n.adapt=5000, n.update=5000, n.iter = 5000, n.chains=nChains) 

stopCluster(cl)
time2run[3] <- round((proc.time()[3]-t.start[3])/(60), 1)
cat(paste("Process time =", time2run[3], "minutes"))

# 1.2 mins for 5000 all
#------------------------------------------------------------------------------
# 4)  Fit MTE - SSUL
#------------------------------------------------------------------------------
cl <- makeCluster(nChains)
t.start[4]<-proc.time()[3]

fitMTE.SSUL <-jags.parfit(cl, data=dat.list, params=c("E", "El", "Eh", "Tl", "Th", "mu0"), model=modelSurvMTE.SSUL, n.adapt=5000, n.update=5000, n.iter = 5000, n.chains=nChains) 

stopCluster(cl)
time2run[4] <- round((proc.time()[3]-t.start[4])/(60), 1)
cat(paste("Process time =", time2run[4], "minutes"))

#------------------------------------------------------------------------------
# 5) Fit MTE - SSU with one location
#------------------------------------------------------------------------------
dat.list_1L<-dat.list
dat.list_1L$nL <- 1
dat.list_1L$loc <- rep(1, length(dat$Dead_before_L3))

cl <- makeCluster(nChains)
t.start[3]<-proc.time()[3]

fitMTE.SSU_1L <-jags.parfit(cl, data=dat.list_1L, params=c("E", "Eh", "Th", "mu0"), model=modelSurvMTE.SSU, n.adapt=5000, n.update=5000, n.iter = 5000, n.chains=nChains) 

stopCluster(cl)
time2run[3] <- round((proc.time()[3]-t.start[3])/(60), 1)
cat(paste("Process time =", time2run[3], "minutes"))

# 1.8 mins
#------------------------------------------------------------------------------
# Save output
#------------------------------------------------------------------------------

# save.image("workspaces/ModelFits_SurvAbnormalAsDead_20190418.RData")
# Note: THese output files were not uploaded to GitHub due to size limitations
# Questions or requests can be directed to the authors.