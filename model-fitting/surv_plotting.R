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
 
# File for looking at the output of model fits, from surv_fitting.R


library(gplots)
library(dclone)

# Function to extract posterior mode
getmode <- function(v) {
	uniqv <- unique(v)
	uniqv[which.max(tabulate(match(v, uniqv)))]
}

###############################################################################
# Separate temperature fits
###############################################################################

# Note that output files were not uploaded to GitHub, but can be generated
# by surv_fitting.R
# load("workspaces/ModelFits_SurvAbnormalAsDead_20190418.RData")

loc <- levels(dat$Location)

#------------------------------------------------------------------------------
# Convergence
#------------------------------------------------------------------------------


Rhat <- gelman.diag(fit)
fit.unlisted <- as.mcmc(rbind(fit[[1]], fit[[2]], fit[[3]], fit[[4]], fit[[5]], fit[[6]], fit[[7]], fit[[8]]))

out <- summary(fit)
o <- colnames(fit[[1]])

# # Option 1: Mean and quantiles of posterior distribution
# mu<-list(
# 	matrix(out[[1]][which(o=="mu[1,1]"):which(o=="mu[10,4]"),1], nrow=dat.list$nT, ncol=dat.list$nL),
# 	matrix(out[[2]][which(o=="mu[1,1]"):which(o=="mu[10,4]"),1], nrow=dat.list$nT, ncol=dat.list$nL),
# 	matrix(out[[2]][which(o=="mu[1,1]"):which(o=="mu[10,4]"),5], nrow=dat.list$nT, ncol=dat.list$nL)
# )

# Option 2: Mode and HPD interval
hpd <- HPDinterval(fit.unlisted)
m <- apply(fit.unlisted, 2, getmode)
mu <- list(
	mode = matrix(m[which(o=="mu[1,1]"):which(o=="mu[10,4]")],nrow = dat.list$nT, ncol = dat.list$nL),
	lower =  matrix(hpd[which(o=="mu[1,1]"):which(o=="mu[10,4]"),'lower'], nrow = dat.list$nT, ncol = dat.list$nL),
	upper =  matrix(hpd[which(o=="mu[1,1]"):which(o=="mu[10,4]"),'upper'], nrow = dat.list$nT, ncol = dat.list$nL)
)

for(k in 1:3){
	rownames(mu[[k]]) <- dat.list$T.obs
	colnames(mu[[k]]) <- loc
	mu[[k]][cbind(c(2,8), rep(which(loc=="Sheep_Mountain"), 2))] <- NA # Set NA to the combinations that had zero data
}

muR <- matrix(Rhat[[1]][which(o=="mu[1,1]"):which(o=="mu[10,4]"),1], nrow=dat.list$nT, ncol=dat.list$nL)

# Create table for plotting/looking at output
write.csv(data.frame(
	Location = rep(c("Colorado", "Alberta", "South Yukon", "North Yukon"), each = dat.list$nT),
	Temperature = rep(dat.list$T.obs, 4),
	Mode = c(mu$mode),
	Lower = c(mu$lower),
	Upper = c(mu$upper),
	Rhat = c(muR)), file = "estimates/muTemp_HPD.csv")

#------------------------------------------------------------------------------
# Plotting temp relationships
#------------------------------------------------------------------------------

# Color points red if they did not converge
colR <- matrix(as.numeric(muR>1.1)+1, nrow=dat.list$nT, ncol=dat.list$nL)

#-----------
# Plot: mortality rate
par(mfrow=c(2,2), mar=c(3,3,2,1), oma=c(2,2,1,0))
for(j in 1:dat.list$nL){
	plotCI(matrix(rep(sort(unique(dat$Temperature)), dat.list$nL), dat.list$nT, dat.list$nL), mu[[1]], li=mu[[2]], ui=mu[[3]], col=grey(0.8), gap=0, main=levels(dat$Location)[j], bty="l", las=1, xlab="", ylab="")
	plotCI(sort(unique(dat$Temperature)), mu[[1]][,j], li=mu[[2]][,j], ui=mu[[3]][,j], add=TRUE, pch=19, cex=0.7, gap=0, col=colR[,j])
}
mtext(side=1, expression(paste("Temperature (", degree, "C)")), outer=TRUE)
mtext(side=2, "Mortality rate (per day)", outer=TRUE)
	
#-----------
# mortality rate (better for comparison)
cols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')
pdf(file="TI_mortRate.pdf", width = 6.5, height = 5, pointsize = 10)
par(mfrow=c(1,1), mar=c(4,4,2,1), oma=c(0,0,0,0))

# par(mfrow=c(1,1), mar=c(4,4,2,1), oma=c(0,0,0,0))
barplot2(t(mu[[1]]), beside=TRUE, plot.ci=TRUE, ci.l = t(mu[[2]]), ci.u = t(mu[[3]]), ylab="Mortality rate (per day)", xlab=expression(paste("Temperature (", degree, "C)")), col=cols, ci.width = 0.2, las=1, ylim=c(0, 0.4))
abline(h = 0)
legend('top', fill=cols, legend=levels(dat$Location), ncol=4, bty="n")

dev.off()
#------------------------------------------------------------------------------
# Table of estimates
#------------------------------------------------------------------------------

mu.table<-data.frame(
	temperature = sort(unique(dat$Temperature)), 
	Colorado = rep(NA, dat.list$nT),
	Sheep_Mountain = rep(NA, dat.list$nT),
	Sheep_River = rep(NA, dat.list$nT),
	Dawson = rep(NA, dat.list$nT)
	)
for (i in 1:dat.list$nT) {
	for (j in 1:dat.list$nL) {
		mu.table[i,j+1] <- paste(formatC(round(mu[[1]][i, j], 4), 4, format="f"), " (", formatC(round(mu[[2]][i, j], 4), 4, format="f"), ", ", formatC(round(mu[[3]][i, j], 4), 4, format="f"), ")", sep="")
	}}
# write.csv(mu.table, file="estimates/TempMuEstimates_nov_abnormAsDead.csv")

NLL <- logLikSurv(params = out[[1]][,1], dat.list = dat.list, mod = "IT")
NLL[[1]]

###############################################################################
# MTE fits
###############################################################################
MTEmod <- "SSU"
# Select which model to examine
# Options: fitMTE.A, fitMTE.SSL, fitMTE.SSUL

loc <- levels(dat$Location)

if(MTEmod == "SSUL"){
	fitMTE <- fitMTE.SSUL
	keepChains <- 1:nChains
	np <- 6 #Number of parameters for the MTE model
	parSym <- c("mu0", "E", "Eh", "El", "Th", "Tl")

	} else if (MTEmod == "SSU"){
	fitMTE <- fitMTE.SSU
	keepChains <- 1:nChains
	np <- 4 # Number of parameters for the MTE model
	parSym <- c("mu0", "E", "Eh", "Th")

} else if (MTEmod == "A"){
	fitMTE <- fitMTE.A
	keepChains <- 1:nChains
	np <- 2 #Number of parameters for the MTE model
	parSym <- c("mu0", "E")

} else if (MTEmod == "SSU_1L"){
	fitMTE <- fitMTE.SSU_1L
	keepChains <- 1:nChains
	np <- 4 #Number of parameters for the MTE model
	parSym <- c("mu0", "E", "Eh", "Th")
	loc <- "All"
}


# #------
# # Keep only chains that converged: (ALL CHAINS CONVERGED; IGNORE)
# fitMTE2 <- list(); length(fitMTE2) <- length(keepChains)
# for(i in 1:length(keepChains)) fitMTE2[[i]] <- fitMTE[[keepChains[i]]]
# fitMTE <- as.mcmc(fitMTE2)

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
# i <- which(oMTE=='mu0[1]')
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
# 
#------------------------------------------------------------------------------
# MTE relationships
#------------------------------------------------------------------------------
T.all<-seq(0, 40, 0.1)

paramsMTE <- list(
	mode = matrix(NA, nrow = np, ncol = dat.list$nL, dimnames=list(parSym, loc)),
	li = matrix(NA, nrow = np, ncol = dat.list$nL, dimnames=list(parSym, loc)),
	ui = matrix(NA, nrow = np, ncol = dat.list$nL, dimnames=list(parSym, loc)))
MTER <- matrix(NA, nrow = np, ncol = dat.list$nL, dimnames=list(parSym, loc))

# if(MTEmod == "SSU_1L"){
# 	for(i in 1:np){
# 		paramsMTE[[1]][i,] <- outMTE[[1]][which(oMTE==parSym[i]),1]
# 		paramsMTE[[2]][i,] <- outMTE[[2]][which(oMTE==parSym[i]),1]
# 		paramsMTE[[3]][i,] <- outMTE[[2]][which(oMTE==parSym[i]),5]
# 		MTER[i,] <- RhatMTE[[1]][which(oMTE==parSym[i])]
# 	}
# } else {

# Using HPD estimates
hpdMTE <- HPDinterval(mcmc(fitMTE.unlisted))
mMTE <- apply(fitMTE.unlisted, 2, getmode)

for(i in 1:np){
	paramsMTE[[1]][i,] <- mMTE[which(oMTE==paste(parSym[i], "[1]", sep="")):which(oMTE==paste(parSym[i], "[4]", sep=""))]
	paramsMTE[[2]][i,] <- hpdMTE[which(oMTE==paste(parSym[i], "[1]", sep="")):which(oMTE==paste(parSym[i], "[4]", sep="")),1]
	paramsMTE[[3]][i,] <- hpdMTE[which(oMTE==paste(parSym[i], "[1]", sep="")):which(oMTE==paste(parSym[i], "[4]", sep="")),2]
	
	MTER[i,] <- RhatMTE[[1]][which(oMTE==paste(parSym[i], "[1]", sep="")):which(oMTE==paste(parSym[i], "[4]", sep=""))]
}
# }

# Table for Supplement
write.csv(data.frame(
	Hyperparam = rep(c("mu0", "E", "Eh","Th"), 4),
	Location = rep(c("Colorado", "Alberta", "South Yukon", "North Yukon"), each =4),
	Mode = c(paramsMTE$mode),
	Lower = c(paramsMTE$li),
	Upper = c(paramsMTE$ui),
	Rhat = c(MTER)), file = "estimates/TableS4_mu.csv")

# Table for export
MTE.table<-data.frame(
	parameter = parSym, 
	Colorado = rep(NA, np),
	Sheep_River = rep(NA, np),
	Sheep_Mountain = rep(NA, np),
	Dawson = rep(NA, np)
)
for (i in 1:np) {
	for (j in 1:dat.list$nL) {
		MTE.table[i,j+1] <- paste(formatC(round(paramsMTE[[1]][i, j], 4), 4, format="f"), " (", formatC(round(paramsMTE[[2]][i, j], 4), 4, format="f"), ", ", formatC(round(paramsMTE[[3]][i, j], 4), 4, format="f"), ")", sep="")
	}}
# write.csv(MTE.table, file="estimates/SSUMuEstimates_abnormalAsDead_Apr18.csv")

#-----------------
# Plot MTE relationships on top of temp dependent tau

par(mfrow=c(2,2), mar=c(3,3,2,1), oma=c(2,2,1,0))
for(j in 1:dat.list$nL){
	plotCI(sort(unique(dat$Temperature)), mu[[1]][,j], li=mu[[2]][,j], ui=mu[[3]][,j], gap=0, bty="l", las=1, xlab="", ylab="", pch=21, col=cols[j], pt.bg="white", lwd=1.2, ylim=c(0, 0.35)) #max(tau[[3]], na.rm=TRUE)
	lines(T.all, predictMTE(Temp = T.all, params = paramsMTE.SSUL[[1]][,j], fn = "SSUL"))
	lines(T.all, predictMTE(Temp = T.all, params = paramsMTE.SSU[[1]][,j], fn = "SSU"), lty=3, lwd=1.5)
	# lines(T.all, predictMTE(Temp = T.all, params = paramsMTE.A[[1]][,j], fn = "A"), lty=2)
	mtext(side=3, levels(dat$Location)[j], line=0.5, col=cols[j], font=2)
}
mtext(side=1, expression(paste("Temperature (", degree, "C)")), outer=TRUE)
mtext(side=2, "Mortality rate (per day)", outer=TRUE)

legend("top", lty = c(2,3,1), lwd=c(1, 1.5, 1), c("Arrhenius", "SS-upper", "SS-upper & lower"), bty="n")

#------------------------------------------------------------------------------
# Parameters from different locations
#------------------------------------------------------------------------------
# Color points red if they did not converge
colRMTE <- matrix(as.numeric(MTER>1.1)+1, nrow=6, ncol=dat.list$nL)

locShort <- c("Col", "ShR", "ShM", "Daw")

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


NLL <- logLikSurv(params = outMTE[[1]][,1], dat.list = dat.list, mod = MTEmod)
NLL[[1]]

#------------------------------------------------------------------------------
# BestMod params table
#------------------------------------------------------------------------------

# tab4R0 <- data.frame(param = rep("mu", length(paramsMTE[[1]])), hyperparam = rep(rownames(paramsMTE[[1]]), 4), location = rep(loc, each = 4), mean = as.numeric(paramsMTE[[1]]), lower = as.numeric(paramsMTE[[2]]), upper = as.numeric(paramsMTE[[3]]))
# write.csv(tab4R0, file = "estimates/bestModParams_abnormalAsDead.csv")
