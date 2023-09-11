###############################################################################
# This code produces figures for the paper 
#
# Local thermal adaptation and local temperature regimes
# drive the performance of parasitic helminths under climate
# change: the case of Marshallagia marshalli from wild
# ungulates 
#
# by Aleuy et al. (2023, Global Change Biology)
#
# Code written by Steph Peacock <stephanie.j.peacock at gmail.com>
# June 26 2023
#
###############################################################################

library(gplots)
library(dclone)

source("functions.R")
# cols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')
cols <- c('#d7191c','#fdae61','#abd9e9','#2c7bb6') # Revised location colors based on Reviewer's comments

################################################################################
# Data (Fig. 4)
################################################################################

dat <- read.csv("data/Marsh_data.csv")

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
#-----------------------------------------------------------------------------
# Development time to L3 
#-----------------------------------------------------------------------------

# Create true location data
dat$Location <- as.character(dat$Location2)
dat$Location[dat$Location == "Sheep_River_ShM" | dat$Location == "Sheep_River_Dawson" | dat$Location == "Sheep_River_Colorado"] <- "Sheep_River"

# Order factor from south to north
dat$Location <- factor(dat$Location, levels=c("Colorado", "Sheep_River", "Sheep_Mountain", "Dawson")) # Sheep_Mountain = South Yukon and Sheep_River = Alberta

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


#-----------------------------------------------------------------------------
# Plot data
#-----------------------------------------------------------------------------
quartz(width = 6.3, height = 3, pointsize = 10)
par(mfrow = c(1, 2), mar = c(4, 2, 1, 1), oma = c(0, 3, 1, 0))

dat1 <- subset(dat, dat$Temperature < 20)
boxplot(dat1$L3_ocurrence ~ dat1$Location + dat1$Temperature, col = paste0(cols, 50), border = cols, at = c(1:4, 6:9, 11:14, 16:19), las = 1, xaxt = "n", xlab = "Temperature", ylab = "")
axis(side = 1, at = c(5, 10, 15), labels = FALSE)
axis(side = 1, at = c(2.5, 7.5, 12.5, 17.5), tck = 0, labels = paste0(sort(unique(dat1$Temperature)), "˚C"))
mtext(side = 2, line = 3.5, "Days to L3")
mtext(side = 3, adj = 0, "a)", line = 0.5)
legend("topright", fill = cols, border = NA, legend = c("Colorado", "Alberta", "South Yukon", "North Yukon"), bty = "n")

dat2 <- subset(dat, dat$Temperature >= 20)
boxplot(dat2$L3_ocurrence ~ dat2$Location + dat2$Temperature, col = paste0(cols, 50), border = cols, at = c(1:4, 6:9, 11:14, 16:19, 21:24), las = 1, xaxt = "n", xlab = "Temperature", ylab = "")
axis(side = 1, at = c(5, 10, 15, 20), labels = FALSE)
axis(side = 1, at = c(2.5, 7.5, 12.5, 17.5, 22.5), tck = 0, labels = paste0(sort(unique(dat2$Temperature))[1:5], "˚C"))
mtext(side = 3, adj = 0, "b)", line = 0.5)

# Extra plot for Sheep River sites
dat3 <- subset(dat, dat$Location == "Sheep_River" & dat$Temperature < 20)
dat3$Location2 <- factor(dat3$Location2, levels= c("Sheep_River_Colorado", "Sheep_River_ShM", "Sheep_River_Dawson"))
boxplot(dat3$L3_ocurrence ~ dat3$Location2 + dat3$Temperature, at = c(1:3, 5:7, 9:11, 13:15), las = 1, xaxt = "n", xlab = "Temperature", ylab = "", col = paste0(cols[c(1,3:4)], 50), border = cols[c(1,3:4)])
axis(side = 1, at = c(4,8,12), labels = FALSE)
axis(side = 1, at = c(2,6,10,14), tck = 0, labels = paste0(sort(unique(dat3$Temperature)), "˚C"))
mtext(side = 3, adj = 0, "a)", line = 0.5)
legend("topright", fill = cols[c(1,3:4)], border = NA, legend = c("Colorado", "South Yukon", "North Yukon"), title = "Alberta comparison to:", bty = "n")

dat4 <- subset(dat, dat$Location == "Sheep_River" & dat$Temperature >= 20)
dat4$Location2 <- factor(dat4$Location2, levels= c("Sheep_River_Colorado", "Sheep_River_ShM", "Sheep_River_Dawson"))
boxplot(dat4$L3_ocurrence ~ dat4$Location2 + dat4$Temperature, col = paste0(cols[c(1,3,4)], 50), border = cols[c(1,3,4)], at = c(1:3, 5:7, 9:11, 13:15, 17:19), las = 1, xaxt = "n", xlab = "Temperature", ylab = "")
axis(side = 1, at = c(4,8,12,16), labels = FALSE)
axis(side = 1, at = c(2,6,10,14,18), tck = 0, labels = paste0(sort(unique(dat2$Temperature))[1:5], "˚C"))
mtext(side = 3, adj = 0, "b)", line = 0.5)

################################################################################
# MTE predictions over temperature (Fig. 5)
###############################################################################

# Vector of dummmy temperatures for plotting
T.all <- seq(-30, 40, 0.1)

#------------------------------------------------------------------------------
# Import best-fit parameter values
#------------------------------------------------------------------------------

paramsRaw <- read.csv("output/estimates/bestMTEParams_abnormalAsDead.csv")

# Create named list for parameters
params <- list(); length(params) <- 4; names(params) <- c("Colorado", "Sheep_River", "Sheep_Mountain", "Dawson") # Order south to north
for(i in 1:4){ # for each location
	params[[i]] <- list(
		tau = paramsRaw$mean[paramsRaw$param=="tau" & paramsRaw$location == names(params)[i]],
		mu = paramsRaw$mean[paramsRaw$param=="mu" & paramsRaw$location == names(params)[i]]
	)
	
	names(params[[i]]$tau) <- as.character(paramsRaw$hyperparam[paramsRaw$param=="tau" & paramsRaw$location == names(params)[i]])
	names(params[[i]]$mu) <- as.character(paramsRaw$hyperparam[paramsRaw$param=="mu" & paramsRaw$location == names(params)[i]])

	} # end location

#------------------------------------------------------------------------------
# Predict tau and mu
#------------------------------------------------------------------------------

# Update Jan 26, 2022: Because of relationships between parameters (in particular Th and Eh) need to calculate MTE curves for each MCMC draw
# Create vector to store predicted demographic parameters from estimated hyperparameters

parPred <- list(
	tau = array(NA, dim = c(4, length(T.all), 3), dimnames = list(c("Colorado", "Alberta", "SouthYukon", "NorthYukon"), T.all, c("mean", "lower", "upper"))),
	mu = array(NA, dim = c(4, length(T.all), 3), dimnames = list(c("Colorado", "Alberta", "SouthYukon", "NorthYukon"), T.all, c("mean", "lower", "upper")))
)
#---
# Tau
#---

tau.mcmc <- readRDS("estimates/MCMC_tauMTE.rds") # Includes only chains that converged (n = 5)

# Note that locations were misnumbered when fitting tau models as c("Colorado", "Sheep_Mountain", "Sheep_River", "Dawson"); need to fix this
locNumTau <- c(1,3,2,4)

o <- colnames(tau.mcmc[[1]])
z <- cbind(sample(1:length(tau.mcmc), size = 1000, replace = TRUE), sample(1:nrow(tau.mcmc[[1]]), size = 1000, replace = TRUE))

for(i in 1:4){
	for(j in 1:length(T.all)){
		tau.ij <- numeric(1000)
		for(k in 1:1000){
			a <- tau.mcmc[[z[k,1]]][cbind(z[k,2], which(o == paste0('tau0[', locNumTau[i], ']')))]
			E <- tau.mcmc[[z[k,1]]][cbind(z[k,2], which(o == paste0('E[', locNumTau[i], ']')))]
			El <- tau.mcmc[[z[k,1]]][cbind(z[k,2], which(o == paste0('El[', locNumTau[i], ']')))]
			Eh <- tau.mcmc[[z[k,1]]][cbind(z[k,2], which(o == paste0('Eh[', locNumTau[i], ']')))]
			Tl <- tau.mcmc[[z[k,1]]][cbind(z[k,2], which(o == paste0('Tl[', locNumTau[i], ']')))]
			Th <- tau.mcmc[[z[k,1]]][cbind(z[k,2], which(o == paste0('Th[', locNumTau[i], ']')))]
			
			tau.ij[k] <- a * exp(1 * E/(8.62*10^-5) * (1/(T.all[j] + 273.15) - 1/(15 + 273.15))) * (1 + exp(El/(8.62*10^-5) * (1/(T.all[j] + 273.15) - 1/(Tl + 273.15))) + exp(Eh/(8.62*10^-5) * (-1/(T.all[j] + 273.15) + 1/(Th + 273.15))))
		}
		
		parPred$tau[i, j, 1] <- mean(tau.ij)
		parPred$tau[i, j, 2:3] <- quantile(tau.ij, c(0.025, 0.975))
	}}

# plot(T.all, parPred$tau[2,,1], "l")

#---
# Mu
#---

mu.mcmc <- readRDS("estimates/MCMC_muMTE.rds")
o <- colnames(mu.mcmc[[1]])
z <- cbind(sample(1:length(mu.mcmc), size = 1000, replace = TRUE), sample(1:nrow(mu.mcmc[[1]]), size = 1000, replace = TRUE))

for(i in 1:4){
	for(j in 1:length(T.all)){
		mu.ij <- numeric(1000)
		for(k in 1:1000){
			a <- mu.mcmc[[z[k,1]]][cbind(z[k,2], which(o == paste0('mu0[', i, ']')))]
			E <- mu.mcmc[[z[k,1]]][cbind(z[k,2], which(o == paste0('E[', i, ']')))]
			Eh <- mu.mcmc[[z[k,1]]][cbind(z[k,2], which(o == paste0('Eh[', i, ']')))]
			Th <- mu.mcmc[[z[k,1]]][cbind(z[k,2], which(o == paste0('Th[', i, ']')))]
		
			mu.ij[k] <- a * exp(-1 * E/(8.62*10^-5) * (1/(T.all[j] + 273.15) - 1/(15 + 273.15))) * (1 + exp(Eh/(8.62*10^-5) * (-1/(T.all[j] + 273.15) + 1/(Th + 273.15))))
		}
		
		parPred$mu[i, j, 1] <- mean(mu.ij)
		parPred$mu[i, j, 2:3] <- quantile(mu.ij, c(0.025, 0.975))
	}}


# parPred <- list(
# 	tau = cbind(
# 		predictMTE(T.all, params = params[[1]]$tau, fn = "SSUL"),
# 		predictMTE(T.all, params = params[[2]]$tau, fn = "SSUL"),
# 		predictMTE(T.all, params = params[[3]]$tau, fn = "SSUL"),
# 		predictMTE(T.all, params = params[[4]]$tau, fn = "SSUL")
# 	),
# 	mu = cbind(
# 		predictMTE(T.all, params = params[[1]]$mu, fn = "SSU"),
# 		predictMTE(T.all, params = params[[2]]$mu, fn = "SSU"),
# 		predictMTE(T.all, params = params[[3]]$mu, fn = "SSU"),
# 		predictMTE(T.all, params = params[[4]]$mu, fn = "SSU")
# 	)
# )
plot(T.all, parPred$mu[2,,1], "l")

# saveRDS(parPred, "estimates/parPred.rds")
# parPred <- readRDS("estimates/parPred.rds")

#------------------------------------------------------------------------------
# Figure 5
#------------------------------------------------------------------------------

quartz(width = 6.3, height = 4.5, pointsize = 10)
par(mfrow = c(2, 2), mar = c(4, 4, 1, 1), oma = c(0, 0, 1, 0))


# Tau parameters
plotCI(1:4, 
			 paramsRaw$mean[paramsRaw$param == "tau" & paramsRaw$hyperparam == "El"], 
			 li = paramsRaw$lower[paramsRaw$param == "tau" & paramsRaw$hyperparam == "El"],
			 ui = paramsRaw$upper[paramsRaw$param == "tau" & paramsRaw$hyperparam == "El"],
			 col= cols, pch = 19, gap = 0,
			 xlim = c(0.5, 9.5), ylim = range(paramsRaw[which(paramsRaw$param == "tau" & paramsRaw$hyperparam %in% c("El", "Eh")), 4:6]), xaxt = "n", las = 1, ylab = "Inactivation energy", xlab = "")

plotCI(6:9, 
			 paramsRaw$mean[paramsRaw$param == "tau" & paramsRaw$hyperparam == "Eh"], 
			 li = paramsRaw$lower[paramsRaw$param == "tau" & paramsRaw$hyperparam == "Eh"],
			 ui = paramsRaw$upper[paramsRaw$param == "tau" & paramsRaw$hyperparam == "Eh"],
			 col= cols, pch = 19, gap = 0,
			 add=TRUE)
abline(v = 5, lty = 2)
axis(side = 1, at = 5, labels = FALSE)
axis(side = 1, at = c(2.5, 7.5), labels= c(expression(paste("Lower (", E[L], ")", sep = "")), expression(paste("Upper (", E[H], ")", sep = ""))), tck = 0)
mtext(side = 3, line = 0.5, "a)", adj = 0)

# Tau over temp
plot(T.all, parPred$tau[1, , 1], "n", ylim = c(0, 80), xlim = c(5,  40), xaxs = "i", las = 1, xlab = expression(paste("Temperature (", degree, "C)", sep = "")), yaxs = "i", ylab = expression(paste("Development time (", tau[L], "; days)")))
for(i in 1:4){
	# polygon(x = c(T.all, rev(T.all)), y = c(parPred$tau[i, , 2], rev(parPred$tau[i, , 3])), border = NA, col = paste0(cols[i], 30))
	lines(T.all, parPred$tau[i, , 1], col = cols[i])
}
mtext(side = 3, line = 0.5, "b)", adj = 0)

# 

# Mu parameters
par(mar = c(4, 13, 1, 1))
plotCI(1:4, 
			 paramsRaw$mean[paramsRaw$param == "mu" & paramsRaw$hyperparam == "Eh"], 
			 li = paramsRaw$lower[paramsRaw$param == "mu" & paramsRaw$hyperparam == "Eh"],
			 ui = paramsRaw$upper[paramsRaw$param == "mu" & paramsRaw$hyperparam == "Eh"],
			 col= cols, pch = 19, gap = 0,
			 xlim = c(0.5, 4.5), ylim = range(paramsRaw[which(paramsRaw$param == "mu" & paramsRaw$hyperparam == "Eh"), 4:6]), xaxt = "n", las = 1, ylab = expression(paste("Upper inactivation energy (", E[H], ")", sep = "")), xlab = "Location")

mtext(side = 3, line = 0.5, "c)", adj = 0)

legend(-5.6, 7.5, fill= cols, border = NA, legend = c("Colorado", "Alberta", "South Yukon", "North Yukon"), bty = "n", xpd = NA)

# Mu over temp
par(mar = c(4, 4, 1, 1))

plot(T.all, parPred$mu[1, , 1], "n", ylim = c(0, 0.3), xaxs = "i", las = 1, xlab = expression(paste("Temperature (", degree, "C)", sep = "")), yaxs = "i", ylab = expression(paste("Mortality rate (", mu[L], "; days", {}^-1, ")")), xlim = c(20,  40))
for(i in 1:4){
	# polygon(x = c(T.all, rev(T.all)), y = c(parPred$mu[i, , 2], rev(parPred$mu[i, , 3])), border = NA, col = paste0(cols[i], 30))
	lines(T.all, parPred$mu[i, , 1], col = cols[i])
}
mtext(side = 3, line = 0.5, "d)", adj = 0)

################################################################################
# Potential performance (Fig. 6)
###############################################################################

dumTemp <- seq(-5, 40, 0.1)
R0pred <- list(); length(R0pred) <- 3
# For three different levels of rhoH
for(j in 1:3){
	R0pred[[j]] <- cbind( # columns = location-specific parameters
		R0calc(temp = dumTemp, params = params[[1]], rhoH = c(0.1, 0.01, 0.001)[j]),
		R0calc(temp = dumTemp, params = params[[2]], rhoH = c(0.1, 0.01, 0.001)[j]),
		R0calc(temp = dumTemp, params = params[[3]], rhoH = c(0.1, 0.01, 0.001)[j]),
		R0calc(temp = dumTemp, params = params[[4]], rhoH = c(0.1, 0.01, 0.001)[j]))
	colnames(R0pred[[j]]) <- dimnames(params) 
	}


# Default plots (main text) are rhoH = 0.01 (i.e., j = 2)
j <- 2
quartz(width = 3.2, height = 2.8, pointsize = 10)
par(mar = c(4, 4.5, 1, 1))

plot(dumTemp, R0pred[[j]][, 1], "n", ylim = c(0, 1), xaxs = "i", las = 1, xlab = expression(paste("Temperature (", degree, "C)", sep = "")), yaxs = "i", ylab = expression(paste("Potential performance (", R[0]/C, ")")))
for(i in 1:4){
	lines(dumTemp,  R0pred[[j]][, i], col = cols[i], lwd = 1.5)
	boop <- which(R0pred[[j]][, i] == max(R0pred[[j]][, i]))
	segments(x0 = dumTemp[boop], x1 = dumTemp[boop], y0 = 0, y1 = R0pred[[j]][boop, i], col = cols[i])
	print(dumTemp[boop])
}
legend("topright", lwd = 1.2, col = cols, legend = c("Colorado", "Alberta", "South Yukon", "North Yukon"), bty = "n", cex = 0.8)

quartz(width = 6.3, height = 2.8, pointsize = 10)
par(mfrow =c (1,2), mar = c(4, 4.5, 1.5, 1))
for(j in c(1,3)){
	plot(dumTemp, R0pred[[j]][, 1], "n", ylim = c(0, 1), xaxs = "i", las = 1, xlab = expression(paste("Temperature (", degree, "C)", sep = "")), yaxs = "i", ylab = expression(paste("Potential performance (", R[0]/C, ")")))
	for(i in 1:4){
		lines(dumTemp,  R0pred[[j]][, i], col = cols[i], lwd = 1.5)
		boop <- which(R0pred[[j]][, i] == max(R0pred[[j]][, i]))
		segments(x0 = dumTemp[boop], x1 = dumTemp[boop], y0 = 0, y1 = R0pred[[j]][boop, i], col = cols[i])
		print(dumTemp[boop])
	}
	if(j == 3) legend("topright", lwd = 1.2, col = cols, legend = c("Colorado", "Alberta", "South Yukon", "North Yukon"), bty = "n", cex = 0.8)
	mtext(side = 3, line = 0.5, adj = 0, c("a)", NA, "b)")[j])
}

################################################################################
# Sesaonal temperature and performance (Fig. 7)
###############################################################################
DOY <- c(1:365)
xDate <- as.Date(paste("2100", DOY, sep = "-"), format = "%Y-%j")

tempLocChange <- readRDS(file = "data/climate-change/tempLocChange.rds")

y_rcp <- readRDS(file = "data/climate-change/tChange.rds")

parChange <- readRDS(file = "data/climate-change/parChange.rds")

R0.CC <- array(NA, dim = c(365, 4, 3), dimnames = list(NULL, LL[1:4], c("current", "rcp26", "rcp85")))

for(i in 1:4){ # For four locations, use parameters from those experiments
	for(j in 1:3){
		R0.CC[, i, j] <- R0calc(temp = tempLocChange[, i, j], params = params[[i]], rhoH = 0.01)
	}}

#------------------------------------------------------------------------------
# Figure
#------------------------------------------------------------------------------
quartz(width = 7.5, height = 5.5, pointsize =10)
par(mfcol = c(3,4), mar = c(2,2,2,1), oma = c(3,3,1,0))
for(i in 1:4){
	
	# Temperature change
	
	plot(xDate, rep(1, 365), "n", ylab = "", las = 1, ylim = c(1, 10), xaxs = "i", xlab = "")
	if(i == 1) mtext(side = 2, expression(paste("Temperature change (", degree, "C)", sep = "")), cex = 0.8, line = 3)
	
	for(j in 1:2){ # for the two climate scenarios
		for(m in 1:12){ # for each month
			polygon(
				x = rep(c(as.Date(paste("2100", m, "01", sep = "-")),as.Date(paste("2100", m, c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[m], sep = "-"))), each = 2), 
				y = y_rcp[i, c(1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1)[m], c('tas25', 'tas75', 'tas75', 'tas25'), j],
				col = c("#CCE5F9", "#F7D4D9")[j], 
				border = NA)
			segments(
				x0 = as.Date(paste("2100", m, "01", sep = "-")), 
				x1 = as.Date(paste("2100", m, c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)[m], sep = "-")), 
				y0 = y_rcp[i, c(1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1)[m], 'tas50', j], 
				y1 = y_rcp[i, c(1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 1)[m], 'tas50', j], 
				col = c("#64B2EC", "#E37E8E")[j], 
				lwd = 2)
			
		} # end m
		
		lines(as.Date(paste("2100", DOY, sep = "-"), format = "%Y-%j"), parChange[i, 'ck', j] + parChange[i, 'dk', j] * cos((DOY - parChange[i, 't0', j]) * 2 * pi / 365), col = c(4,2)[j], lwd = 1.5)
		
	} # end j
	mtext(side = 3, line = 1.5, LL[i], cex = 0.8)
	if(i == 1) mtext(side = 3, line = 0.5, adj = 0, "a)", cex = 0.8)
	if(i == 1) legend("topleft", col = c(4,2), legend = c("Low emissions", "High emissions"), lwd = 1.5, bty = "n")
	
	# Temperature regimes
	plot(xDate, tempLocChange[, i, 1], "l", xlab = "", ylab= "", las = 1, ylim = c(-25, 25), xaxs = "i")
	abline(h = 0)
	for(j in 1:2){
		lines(as.Date(paste("2100", DOY, sep = "-"), format = "%Y-%j"), tempLocChange[, i, j + 1], col = c(4,2)[j], lwd = 1.5)
	}
	if(i == 1) mtext(side = 2, expression(paste("Temperature (", degree, "C)", sep = "")), cex = 0.8, line = 3)
	if(i == 1) mtext(side = 3, line = 0.5, adj = 0, "b)", cex = 0.8)
	
	# R0
	plot(xDate, R0.CC[, i, 1], "l", xlab = "", ylab= "", las = 1, ylim = range(R0.CC), xaxs = "i")
	abline(h = 0)
	for(j in 1:2){
		lines(as.Date(paste("2100", DOY, sep = "-"), format = "%Y-%j"), R0.CC[, i, j + 1], col = c(4,2)[j], lwd = 1.5)
	}
	if(i == 1) mtext(side = 2, expression(R[0]/C), cex = 0.8, line = 3)
	if(i == 1) mtext(side = 3, line = 0.5, adj = 0, "c)", cex = 0.8)
	
}

###############################################################################
# Fig. 7
###############################################################################

# see potentialPerformanceCalc.R