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


# File for looking at the output of model fits, from fitting.R

library(gplots)
library(dclone)

# Function to extract posterior mode
getmode <- function(x) {
	d <- density(x)
	d$x[which.max(d$y)]
}

###############################################################################
# Separate temperature fits
###############################################################################

# Load output from developTime_fitting.R
# Note that output files were too large to share via GitHub; if there are issues 
# generating output files or original output is needed email authors
# load("workspaces/AllModelFits_20190123.RData")

loc <- levels(dat$Location)

fit.unlisted <- mcmc(rbind(fit[[1]], fit[[2]], fit[[3]]))

#------------------------------------------------------------------------------
# Convergence
#------------------------------------------------------------------------------

Rhat <- gelman.diag(fit)

out <- summary(fit)
o <- colnames(fit[[1]])

# # Option 1: Mean and quantiles of posterior distribution
# tau <- list(
# 	mean = matrix(out[[1]][which(o=="tau[1,1]"):which(o=="tau[10,4]"),1], nrow = dat.list$nT, ncol = dat.list$nL),
# 	HPD = 
# 	lower = matrix(out[[2]][which(o=="tau[1,1]"):which(o=="tau[10,4]"),1], nrow=dat.list$nT, ncol=dat.list$nL),
# 	upper = matrix(out[[2]][which(o=="tau[1,1]"):which(o=="tau[10,4]"),5], nrow=dat.list$nT, ncol=dat.list$nL)
# )

# Option 2: Mode and HPD interval
hpd <- HPDinterval(fit.unlisted)
m <- apply(fit.unlisted, 2, getmode)
tau <- list(
	mode = matrix(m[which(o=="tau[1,1]"):which(o=="tau[10,4]")],nrow = dat.list$nT, ncol = dat.list$nL),
	lower =  matrix(hpd[which(o=="tau[1,1]"):which(o=="tau[10,4]"),'lower'], nrow = dat.list$nT, ncol = dat.list$nL),
	upper =  matrix(hpd[which(o=="tau[1,1]"):which(o=="tau[10,4]"),'upper'], nrow = dat.list$nT, ncol = dat.list$nL)
)
# saveRDS(tau, "tauTempSpecificEstimates_HPD.rds")


for(k in 1:3){
	rownames(tau[[k]]) <- dat.list$T.obs
	colnames(tau[[k]]) <- loc
	tau[[k]][cbind(c(2,8), rep(which(loc=="Sheep_Mountain"), 2))] <- NA # Set NA to the combinations that had zero data
}

tauR <- matrix(Rhat[[1]][which(o=="tau[1,1]"):which(o=="tau[10,4]"),1], nrow=dat.list$nT, ncol=dat.list$nL)

# Create table for plotting/looking at output
write.csv(data.frame(
	Location = rep(c("Colorado", "South Yukon", "Alberta", "North Yukon"), each = dat.list$nT), #Note that I had the locations South Yukon and Alberta numbered wrong (not S->N) in the data list
	Temperature = rep(dat.list$T.obs, 4),
	Mode = c(tau$mode),
	Lower = c(tau$lower),
	Upper = c(tau$upper),
	Rhat = c(tauR)), file = "estimates/tauTemp_HPD.csv")

# Seems like the low temperatures had convergence issues...look at trace plots
i <- which(o=='sigma_v')
y <- fit[[1]][,i]; for(j in 2:nChains) y <- c(y, fit[[j]][,i])
par(mfrow=c(1,2), mar=c(4,4,2,1))
plot(c(1:dim(fit[[1]])[1]), fit[[1]][,i], type= "n", ylab=o[i], xlab="", ylim=range(y), main="MCMC trace")
for(j in 1:length(fit)) lines(c(1:dim(fit[[j]])[1]), fit[[j]][,i], lty=3, col=j)

plot(density(fit[[1]][,i]), xlim=range(y), main="Density")
for(j in 1:length(fit)) lines(density(fit[[j]][,i]), col=j)
legend('topleft', lwd=1, col=1:nChains, legend=1:nChains, ncol=1, bty="n")
# add prior:
lines(seq(0, 200, 0.1), dexp(seq(0, 200, 0.1), rate=1/100), col=grey(0.8), lwd=2, lty=3)

#------------------------------------------------------------------------------
# Plotting temp relationships
#------------------------------------------------------------------------------

# Color points red if they did not converge
colR <- matrix(as.numeric(tauR>1.1)+1, nrow=dat.list$nT, ncol=dat.list$nL)

#-----------
# Plot: development time
par(mfrow=c(2,2), mar=c(3,3,2,1), oma=c(2,2,1,0))
for(j in 1:dat.list$nL){
	plotCI(matrix(rep(sort(unique(dat$Temperature)), dat.list$nL), dat.list$nT, dat.list$nL), tau[[1]], li=tau[[2]], ui=tau[[3]], col=grey(0.8), gap=0, main=levels(dat$Location)[j], bty="l", las=1, xlab="", ylab="", ylim=c(0, 250))
	plotCI(sort(unique(dat$Temperature)), tau[[1]][,j], li=tau[[2]][,j], ui=tau[[3]][,j], add=TRUE, pch=19, cex=0.7, gap=0, col=colR[,j])
}
mtext(side=1, expression(paste("Temperature (", degree, "C)")), outer=TRUE)
mtext(side=2, "Time to L3 (days)", outer=TRUE)
	
#-----------
# Plot: development rate
par(mfrow=c(2,2), mar=c(3,3,2,1), oma=c(2,2,1,0))
for(j in 1:dat.list$nL){
	plotCI(matrix(rep(sort(unique(dat$Temperature)), dat.list$nL), dat.list$nT, dat.list$nL), 1/tau[[1]], li=1/tau[[3]], ui=1/tau[[2]], col=grey(0.8), gap=0, main=levels(dat$Location)[j], bty="l", las=1, xlab="", ylab="", ylim=c(0, 0.2))
	plotCI(sort(unique(dat$Temperature)), 1/tau[[1]][,j], li=1/tau[[3]][,j], ui=1/tau[[2]][,j], add=TRUE, pch=19, cex=0.7, gap=0, col=colR[,j])
}
mtext(side=1, expression(paste("Temperature (", degree, "C)")), outer=TRUE)
mtext(side=2, expression(paste("Development rate (day", {}^-1, ")")), outer=TRUE)

#-----------
# Development time (better for comparison)
cols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')
pdf(file="TI_devTime.pdf", width = 6.5, height = 5, pointsize = 10)
par(mfrow=c(2,1), mar=c(4,4,2,1), oma=c(0,0,0,0))

# par(mfrow=c(1,1), mar=c(4,4,2,1), oma=c(0,0,0,0))
barplot2(t(tau[[1]]), beside=TRUE, plot.ci=TRUE, ci.l = t(tau[[2]]), ci.u = t(tau[[3]]), ylab="Time to L3 (days)", xlab=expression(paste("Temperature (", degree, "C)")), col=cols, ci.width = 0.2, las=1, ylim=c(0, 300))
abline(h = 0)
legend('top', fill=cols, legend=levels(dat$Location), ncol=4, bty="n")

# par(mfrow=c(1,1), mar=c(4,4,2,1), oma=c(0,0,0,0))
barplot2(log(t(tau[[1]])), beside=TRUE, plot.ci=TRUE, ci.l = log(t(tau[[2]])), ci.u = log(t(tau[[3]])), ylab="log Time to L3 (days)", xlab=expression(paste("Temperature (", degree, "C)")), col=cols, ci.width = 0.2, las=1)
abline(h = 0)
# legend('top', fill=2:5, legend=levels(dat$Location), ncol=4, bty="n")
dev.off()

#------------------------------------------------------------------------------
# Table of estimates
#------------------------------------------------------------------------------

tau.table<-data.frame(
	temperature = sort(unique(dat$Temperature)), 
	Colorado = rep(NA, dat.list$nT),
	Alberta = rep(NA, dat.list$nT),
	SouthYukon = rep(NA, dat.list$nT),
	NorthYukon = rep(NA, dat.list$nT)
	)
for (i in 1:dat.list$nT) {
	for (j in 1:dat.list$nL) {
		tau.table[i,j+1] <- paste(formatC(round(tau[[1]][i, j], 1), 1, format="f"), " (", formatC(round(tau[[2]][i, j], 1), 1, format="f"), ", ", formatC(round(tau[[3]][i, j], 1), 1, format="f"), ")", sep="")
	}}
# write.csv(tau.table, file="TempTauEstimates_HPD.csv")

NLL <- logLikDevelop(params = out[[1]][,1], dat.list = dat.list, mod = "IT")


###############################################################################
# MTE fits
###############################################################################

load("workspaces/ModelFits_A_SSUL_20190124.RData")

MTEmod <- "SSUL"
# Select which model to examine
# Options: fitMTE.A, fitMTE.SSL, fitMTE.SSUL
if(MTEmod == "SSUL"){
	fitMTE <- fitMTE.SSUL
	# keepChains <- c(1, 3:nChains) #January 22 Great! Everything converged if we do that.
	# keepChains <- c(1,3,5) #January 23, revised prior on tau0. Needs longer?
	keepChains <- c(1,2,4,7,8) #January 24 longer run. Chians 3, 5, & 6 had high sigma_v
	np <- 6 #Number of parameters for the MTE model
	parSym <- c("tau0", "E", "Eh", "El", "Th", "Tl")
} else if (MTEmod == "SSL"){
	fitMTE <- fitMTE.SSL
	keepChains <- c(3:6) #January 23 Great!
	np <- 4 #Number of parameters for the MTE model
	parSym <- c("tau0", "E", "El", "Tl")

} else if (MTEmod == "A"){
	fitMTE <- fitMTE.A
	# keepChains <- c(1,2,4,6) #January 23 
	keepChains <- c(1,2,6,8)
	np <- 2 #Number of parameters for the MTE model
	parSym <- c("tau0", "E")

} else if (MTEmod == "SSUL_1L"){
	fitMTE <- fitMTE.SSUL_1L
	keepChains <- c(1,2,4,6,7,8) #January 28 
	np <- 6 #Number of parameters for the MTE model
	parSym <- c("tau0", "E", "Eh", "El", "Th", "Tl")
	loc <- "All"
}


#------
# Keep only chains that converged:
fitMTE2 <- list(); length(fitMTE2) <- length(keepChains)
for(i in 1:length(keepChains)) fitMTE2[[i]] <- fitMTE[[keepChains[i]]]
fitMTE <- mcmc.list(fitMTE2)

#------------------------------------------------------------------------------
# Convergence
#------------------------------------------------------------------------------

RhatMTE <- gelman.diag(fitMTE)

outMTE <- summary(fitMTE)
oMTE <- colnames(fitMTE[[1]])

# Unlist all chains to get total MCMC samples
fitMTE.unlisted <- fitMTE[[1]]
for(j in 2:length(fitMTE)) fitMTE.unlisted <- rbind(fitMTE.unlisted, fitMTE[[j]])

# # Seems like the low temperatures had convergence issues...look at trace plots
# oMTE[which(RhatMTE[[1]][,1] > 1.1)]
# i <- which(oMTE=='tau0')
# par(mfrow=c(1,2), mar=c(4,4,2,1))
# plot(c(1:dim(fitMTE[[1]])[1]), fitMTE[[1]][,i], type= "n", ylab=o[i], xlab="", ylim=range(fitMTE2[,i]), main="MCMC trace")
# for(j in 1:length(fitMTE)) lines(c(1:dim(fitMTE[[j]])[1]), fitMTE[[j]][,i], lty=3, col=keepChains[j])
# 
# plot(density(fitMTE[[1]][,i]), xlim=range(fitMTE2[,i]), main="Density")
# for(j in 1:length(fitMTE)) lines(density(fitMTE[[j]][,i]), col=keepChains[j])
# legend('topleft', lwd=1, col=keepChains, legend=keepChains, ncol=1, bty="n")
# mtext(side=3, outer=TRUE, round(RhatMTE[[1]][i,1], 2))
# 
# 
# # Add prior for tau0
# priorX <- seq(min(fitMTE2[,i]), max(fitMTE2[,i]), length.out=100)
# priorD <- dlnorm(priorX, log(0.1), 1)
# lines(priorX, priorD, lwd=2, col="#00000040", lty=3)

#------------------------------------------------------------------------------
# MTE relationships
#------------------------------------------------------------------------------
T.all<-seq(0, 40, 0.1)
# paramsMTE <- list(
# 	mean = matrix(NA, nrow = np, ncol = dat.list$nL, dimnames=list(parSym, loc)),
# 	li = matrix(NA, nrow = np, ncol = dat.list$nL, dimnames=list(parSym, loc)),
# 	ui = matrix(NA, nrow = np, ncol = dat.list$nL, dimnames=list(parSym, loc)))

paramsMTE <- list(
	mode = matrix(NA, nrow = np, ncol = dat.list$nL, dimnames=list(parSym, loc)),
	li = matrix(NA, nrow = np, ncol = dat.list$nL, dimnames=list(parSym, loc)),
	ui = matrix(NA, nrow = np, ncol = dat.list$nL, dimnames=list(parSym, loc)))

MTER <- matrix(NA, nrow = np, ncol = dat.list$nL, dimnames=list(parSym, loc))

# if(MTEmod == "SSUL_1L"){
# 	for(i in 1:np){
# 		paramsMTE[[1]][i,] <- outMTE[[1]][which(oMTE==parSym[i]),1]
# 		paramsMTE[[2]][i,] <- outMTE[[2]][which(oMTE==parSym[i]),1]
# 		paramsMTE[[3]][i,] <- outMTE[[2]][which(oMTE==parSym[i]),5]
# 		MTER[i,] <- RhatMTE[[1]][which(oMTE==parSym[i])]
# 	}
# } else {
# for(i in 1:np){
# 	paramsMTE[[1]][i,] <- outMTE[[1]][which(oMTE==paste(parSym[i], "[1]", sep="")):which(oMTE==paste(parSym[i], "[4]", sep="")),1]
# 	paramsMTE[[2]][i,] <- outMTE[[2]][which(oMTE==paste(parSym[i], "[1]", sep="")):which(oMTE==paste(parSym[i], "[4]", sep="")),1]
# 	paramsMTE[[3]][i,] <- outMTE[[2]][which(oMTE==paste(parSym[i], "[1]", sep="")):which(oMTE==paste(parSym[i], "[4]", sep="")),5]
# 	
# 	MTER[i,] <- RhatMTE[[1]][which(oMTE==paste(parSym[i], "[1]", sep="")):which(oMTE==paste(parSym[i], "[4]", sep=""))]
# }
# }

# Using HPD estimates
hpdMTE <- HPDinterval(mcmc(fitMTE.unlisted))
mMTE <- apply(fitMTE.unlisted, 2, getmode)

for(i in 1:np){
	paramsMTE[[1]][i,] <- mMTE[which(oMTE==paste(parSym[i], "[1]", sep="")):which(oMTE==paste(parSym[i], "[4]", sep=""))]
	paramsMTE[[2]][i,] <- hpdMTE[which(oMTE==paste(parSym[i], "[1]", sep="")):which(oMTE==paste(parSym[i], "[4]", sep="")),1]
	paramsMTE[[3]][i,] <- hpdMTE[which(oMTE==paste(parSym[i], "[1]", sep="")):which(oMTE==paste(parSym[i], "[4]", sep="")),2]
	
	
	MTER[i,] <- RhatMTE[[1]][which(oMTE==paste(parSym[i], "[1]", sep="")):which(oMTE==paste(parSym[i], "[4]", sep=""))]
}


# Table for export
MTE.table<-data.frame(
	parameter = parSym, 
	Colorado = rep(NA, np),
	Sheep_Mountain = rep(NA, np), # Note reversal of location N-> S
	Sheep_River = rep(NA, np),
	Dawson = rep(NA, np)
)
for (i in 1:np) {
	for (j in 1:dat.list$nL) {
		MTE.table[i,j+1] <- paste(formatC(round(paramsMTE[[1]][i, j], 2), 2, format="f"), " (", formatC(round(paramsMTE[[2]][i, j], 2), 2, format="f"), ", ", formatC(round(paramsMTE[[3]][i, j], 2), 2, format="f"), ")", sep="")
	}}
# write.csv(MTE.table, file="ATauEstimates_Jan24.csv")

# Table for Supplement
write.csv(data.frame(
	Hyperparam = rep(c("tau0", "E", "Eh", "El", "Th", "Tl"), 4),
	Location = rep(c("Colorado","South Yukon", "Alberta",  "North Yukon"), each = 6), # Note location mis-ordered - as done in dat.list
	Mode = c(paramsMTE$mode),
	Lower = c(paramsMTE$li),
	Upper = c(paramsMTE$ui),
	Rhat = c(MTER)), file = "estimates/tauMTE_HPD2.csv")

#-----------------
# Plot MTE relationships on top of temp dependent tau

par(mfrow=c(2,2), mar=c(3,3,2,1), oma=c(2,2,1,0))
for(j in 1:dat.list$nL){
	plotCI(sort(unique(dat$Temperature)), tau[[1]][,j], li=tau[[2]][,j], ui=tau[[3]][,j], gap=0, bty="l", las=1, xlab="", ylab="", ylim=c(0, 50), pch=21, col=cols[j], pt.bg="white", lwd=1.2) #max(tau[[3]], na.rm=TRUE)
	lines(T.all, predictMTE(Temp = T.all, params = paramsMTE.SSUL[[1]][,j], fn = "SSUL"))
	lines(T.all, predictMTE(Temp = T.all, params = paramsMTE.SSL[[1]][,j], fn = "SSL"), lty=3, lwd=1.5)
	lines(T.all, predictMTE(Temp = T.all, params = paramsMTE.A[[1]][,j], fn = "A"), lty=2)
	mtext(side=3, levels(dat$Location)[j], line=0.5, col=cols[j], font=2)
}
mtext(side=1, expression(paste("Temperature (", degree, "C)")), outer=TRUE)
mtext(side=2, "Time to L3 (days)", outer=TRUE)

legend("top", lty = c(2,3,1), lwd=c(1, 1.5, 1), c("Arrhenius", "SS-lower", "SS-upper & lower"), bty="n")

#------------------------------------------------------------------------------
# Parameters from different locations
#------------------------------------------------------------------------------
# Color points red if they did not converge
colRMTE <- matrix(as.numeric(MTER>1.1)+1, nrow=6, ncol=dat.list$nL)

locShort <- c("Col", "ShM", "ShR", "Daw")

par.names <- c(
	expression(paste("Base rate (", tau[0], ")", sep="")), 
	paste("Activation energy (E)"), expression(paste("Upper inact. energy (", E[h], ")", sep="")), 
	expression(paste("Lower inact. energy (", E[l], ")", sep="")), 
	expression(paste("Upper thermal bound (", T[h], ")", sep="")), 
	expression(paste("Lower thermal bound (", T[l], ")", sep="")))

if(MTEmod=="SSUL"){
	par(mfrow=c(2,3), mar=c(3,4.5,2,1), oma=c(2,0,1,0))
} else if(MTEmod=="SSL"){
	par(mfrow=c(2,2), mar=c(3,4,2,1), oma=c(2,0,1,0))
	par.names <- par.names[c(1,2,4,6)]
} else if(MTEmod=="A"){
	par(mfrow=c(1,2), mar=c(3,4,2,1), oma=c(2,0,1,0))
	par.names <- par.names[c(1,2)]}
for(i in 1:np){
	plotCI(1:dat.list$nL, paramsMTE[[1]][i,], li = paramsMTE[[2]][i,], ui = paramsMTE[[3]][i,], xlab="", xaxt="n", ylab=par.names[i], bty="l", las=1, gap=0, pt.bg="white", xlim=c(0.5,4.5), pch = c(19, 21)[colRMTE[i,]], cex=1.5, lwd=1.2, cex.axis=1.2, cex.lab=1.2, xpd=NA, col=cols)
	axis(side=1, at=c(1:dat.list$nL), labels=locShort) 
	if(i==2) abline(h=0.65, lty=3)
	if(i==3|i==4) abline(h=3.25, lty=3)
	mtext(side=3, adj=0, line=0.5, paste(letters[i], ")", sep=""))
}
mtext(side=1, outer=TRUE, "Location", line=0.5)


NLL <- logLikDevelop(params = outMTE[[1]][,1], dat.list = dat.list, mod = "SSUL_1L")
fitMTE$BUGSoutput$DIC

###############################################################################
# Plot fits to data
###############################################################################
cols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')

params = out[[1]][,1]
dat.list = dat.list
mod = "IT"

locations <- c("Colorado","South Yukon", "Alberta",  "North Yukon") # Note location mis-ordered - as done in dat.list
temps <-  c(5, 7, 10, 15, 20, 25, 30, 33, 35, 38)

# For Colorado, plot development time over temperature
i <- 1
ind <- which(dat.list$loc == which(locations == "Colorado"))
length(ind)

x <- jitter(temps[dat.list$temp[ind]], factor = 1.5)

plot(x, dat.list$time[ind], "n", xlab = "Temperature (˚C)", ylab = "Time to L3 (days)")
segments(x0 = x[which(dat.list$censored[ind] == 0)],
				 x1 = x[which(dat.list$censored[ind] == 0)],
				 y0 = dat.list$time[ind[which(dat.list$censored[ind] == 0)]],
				 y1 = dat.list$time_back[ind[which(dat.list$censored[ind] == 0)]],
				 col = paste0(cols[1], 40))

# From data compi;lation:
# To avoid errors, set NA time_back values to 1 
# Note: these do not come into model fitting as they are censored.

points(x, dat.list$time[ind],  pch = 19, col = paste0(cols[1], 40), cex = 0.6)
points(x, dat.list$time_back[ind],  pch = 19, col = paste0(cols[i], 40), cex = 0.6)

# Add IT model prediction
plotCI(temps, tau$mode[, i], li = tau$lower[, i], ui  = tau$upper[, i], pch = 1, add = TRUE, gap = 0)




dat.list$time[which(dat.list$loc == which(locations == "Colorado") & dat.list$temp == 1)]
