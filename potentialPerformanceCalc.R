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
# Code written by Steph Peacock <stephanie.j.peacock at gmail.com>
# April 2023
#
###############################################################################

library(dclone)
library(gplots)

# Calculate R0 curves under three scenarios:
# (1) Same annual temperature curve, same projected change
# (2) Local annual temperature curve, same projected change
# (3) Local annual temperature curve, local projected change

# Locations
locDat <- read.csv("data/Locations.csv")
LL <- c("Colorado", "SheepRiver", "SheepMountain", "Dawson", "Sierra", "Desert")

###############################################################################
# Functions for R0
###############################################################################

source("functions.R")
cols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')
colP <- c(1:6)

DOY <- c(1:365)
xDate <- as.POSIXct(paste(rep(2000, 365), DOY, rep(0, length(DOY)), sep="-"), format="%Y-%j")

###############################################################################
# Load temperature regimes and changes
###############################################################################

tempLoc <- readRDS("data/climate-change/AnnualTempCurves_20210528.rds")

tempLocChange <- readRDS(file = "data/climate-change/tempLocChange.rds")


dim(tempLoc)
dim(tempLocChange)

i <- 2
plot(DOY, tempLoc[, i, 3], "l", ylim = c(-10, 30), main = LL[i])
lines(DOY, tempLocChange[, i, 1], lty = 2, lwd = 2)
lines(DOY, tempLocChange[, i, 2], lty = 2, lwd = 2, col = 4)
lines(DOY, tempLocChange[, i, 3], lty = 2, lwd = 2, col = 2)
lines(DOY, tempLoc[, 7, 3], lty = 3)

# Illustrate scenarios
quartz(width = 3.5, height = 6)
par(mfrow = c(3,1), mar = c(2,5,2,1), oma = c(0,0,0,0))

# a) Common temp, common change
plot(xDate, tempLoc[, 7, 3], "l", las = 1, ylab = "", ylim = range(tempLocChange))
lines(xDate, tempLoc[, 7, 3] + 5, lty = 2)
mtext(side = 3, adj= 0, "A) Common temp, common change", cex = 0.8)

# b) Local temp, common change
plot(xDate, tempLoc[, 7, 3], "n", las = 1, ylab = "", ylim = range(tempLocChange))
for(i in c(1,4)){
	lines(xDate, tempLoc[, i, 3], col = cols[i])
	lines(xDate, tempLoc[, i, 3] + 5, lty = 2, col =cols[i])
}
mtext(side = 3, adj= 0, "B) Local temp, common change", cex = 0.8)
# c) Local temp, local change
plot(xDate, tempLoc[, 7, 3], "n", las = 1, ylab = "", ylim = range(tempLocChange))
for(i in c(1,4)){
	lines(xDate, tempLocChange[, i, 1], col = cols[i])
	lines(xDate, tempLocChange[, i, 3], lty = 2, col =cols[i])
}
mtext(side = 3, adj= 0, "C) Local temp, local change", cex = 0.8)
mtext(side = 2, outer= TRUE, expression(paste("Temperature (", degree, "C)")), line = -2)
################################################################################
# Load parasite MTE parameters
################################################################################

paramsRaw <- read.csv("estimates/bestModParams_abnormalAsDead.csv")
params <- list(); length(params) <- 4; names(params) <- c("Colorado", "Sheep_River", "Sheep_Mountain", "Dawson") # Order south to north
for(i in 1:4){ # for each location
	params[[i]] <- list(
		tau = paramsRaw$mean[paramsRaw$param=="tau" & paramsRaw$location == names(params)[i]],
		mu = paramsRaw$mean[paramsRaw$param=="mu" & paramsRaw$location == names(params)[i]]
	)
	
	names(params[[i]]$tau) <- as.character(paramsRaw$hyperparam[paramsRaw$param=="tau" & paramsRaw$location == names(params)[i]])
	names(params[[i]]$mu) <- as.character(paramsRaw$hyperparam[paramsRaw$param=="mu" & paramsRaw$location == names(params)[i]])
} # end location


################################################################################
# Potential performance projections under scenarios
################################################################################

# Report DL or R0/C?
DL.metric <- FALSE

# Create arrays to store output
R0.all <- array(NA, 
								dim = c(3, 365, 4, 3), 
								dimnames = list(c("scenario1", "scenario2", "scneario3"), NULL, LL[1:4], c("current", "rcp2.6", 'rcp8.5')))

AUC.all <- array(NA, 
								 dim = c(3, 4, 3), 
								 dimnames = list(c("scenario1", "scenario2", "scneario3"), LL[1:4], c("current", "rcp2.6", 'rcp8.5')))

#------------------------------------------------------------------------------
# Scenario 1: Sample annual curve, same projected change
#------------------------------------------------------------------------------

# Projected change = +2 and +5 degrees all year
for(i in 1:4){ # For four locations, use parameters from those experiments
	for(j in 1:3){ # for different climate change scenarios
			R0.all[1, , i, j] <- R0calc(temp = tempLoc[, 7, 3] + c(0, 2, 5)[j], params = params[[i]], DL.out = DL.metric, rhoH = 0.01)
			
			AUC.all[1, i, j] <- sum(R0.all[1, , i, j])
		}}

# i <- 1
# plot(DOY, R0.all[1, , i, 1], "l", main = LL[i])
# lines(DOY, R0.all[1, , i, 2], col = 4)
# lines(DOY, R0.all[1, , i, 3], col = 2)

#------------------------------------------------------------------------------
# Scenario 2: Local annual curve, same projected change
#------------------------------------------------------------------------------

for(i in 1:4){ # For four locations, use parameters from those experiments
	for(j in 1:3){ # for different climate change scenarios
		R0.all[2, , i, j] <- R0calc(temp = tempLoc[, i, 3] + c(0, 2, 5)[j], params = params[[i]], DL.out = DL.metric, rhoH = 0.01)
		
		AUC.all[2, i, j] <- sum(R0.all[2, , i, j])
	}}

#------------------------------------------------------------------------------
# Scenario 3: Local annual curve, local projected change
#------------------------------------------------------------------------------

for(i in 1:4){ # For four locations, use parameters from those experiments
	for(j in 1:3){ # for different climate change scenarios
		R0.all[3, , i, j] <- R0calc(temp = tempLocChange[, i, j], params = params[[i]], DL.out = DL.metric, rhoH = 0.01)
		
		AUC.all[3, i, j] <- sum(R0.all[3, , i, j])
	}}

################################################################################
# Plot
################################################################################
scenName <- c("Common temp., common change", "Local temp., common change", "Local temp., local change")

quartz()
par(mfrow = c(2, 3))
for(s in 1:3){
	plot(1:4-0.1, AUC.all[s, , 1], xlim = c(0.7, 4.3), col = cols, lwd = 2, cex = 1.5, ylim = range(AUC.all), main = paste("Scenario", s), xaxt = "n", xlab = "", las = 1, ylab = expression(paste(sum(R[0]/C))))
	points(1:4, AUC.all[s, , 2], col = cols, lwd = 2, cex = 1.5, pch = 2)
	points(1:4+0.1, AUC.all[s, , 3], col = cols,  cex = 1.5, pch = 17)
	abline(v = seq(1.5, 3.5, 1))
	axis(side = 1, at = c(1:4), LL[1:4])
}

# As a percent
quartz(width = 4, height = 6.5)
par(mfrow = c(3,1), mar = c(1.5,3,1,1), oma = c(6,2,1,7))
for(s in 1:3){
	plot(1:4-0.1, (AUC.all[s, , 1] - AUC.all[s, , 1])/AUC.all[s, , 1], xlim = c(0.7, 4.3), col = cols, lwd = 2, cex = 1.5, xaxt = "n", xlab = "", las = 1, ylab = "", ylim = c(0, 0.25))
	mtext(side = 3, adj = 0, paste0(LETTERS[s], ") ", scenName[s]), cex = 0.8)
	abline(h = 0, lty = 3)
	points(1:4, (AUC.all[s, , 2] - AUC.all[s, , 1])/AUC.all[s, , 1], col = cols, lwd = 2, cex = 1.5, pch = 2)
	points(1:4+0.1, (AUC.all[s, , 3] - AUC.all[s, , 1])/AUC.all[s, , 1], col = cols,  cex = 1.5, pch = 17)
	abline(v = seq(1.5, 3.5, 1))
	
	if(s == 1) legend(4.6, 0.25, pch = c(1, 2, 17), pt.lwd = 2, col = 1, legend = c("Current", "Low emissions", "High emissions"), xpd = NA, pt.cex = 1.5, bty = "n")
}
axis(side = 1, at = c(1:4),, labels = FALSE)
mtext(side = 2, expression(paste("Percent change in ", sum(R[0]/C), " with climate change")), outer = TRUE, line = 0)
text(1:4, -0.03, srt = 90, c("Colorado", "Sheep\nRiver", "Sheep\nMountain", "Dawson"), col = cols, xpd = NA, adj = 1)
mtext(side = 1, line = 5, "Source of larvae")

# Single panel, under high emissions
AUC.all.percent <- AUC.all
for(s in 1:3){
	for(j in 1:3){
		AUC.all.percent[s, , j] <- (AUC.all[s, , j] - AUC.all[s, , 1])/AUC.all[s, , 1]
}}

j <- 3
# quartz(width = 6.3, height = 2.5, pointsize = 10)
par(mar = c(4,6,2,15))
plot(1:4-0.1, AUC.all.percent[1, , j], ylim = extendrange(AUC.all.percent[, , j], f = 0.1), xlim = c(0.7, 4.3), col = cols, lwd = 2, xaxt = "n", xlab = "", las = 1, ylab = "")
polygon(x = c(0,0,5,5), y = c(range(AUC.all.percent[3, , j]), rev(range(AUC.all.percent[3, , j]))), col = "#00000010", border =NA)
points(1:4, AUC.all.percent[2, , j], col = cols, pch = 19)
points(1:4 + 0.1, AUC.all.percent[3, , j], col = cols, pch = 8, lwd = 1.5)
mtext(side = 2, "Percent change in AUC\n under high emissions scenario", line = 3)

u <- par('usr')
abline(v = seq(1.5, 3.5, 1))
axis(side = 1, at = c(1:4),, labels = FALSE)
text(1:4, u[3] - 0.2*(u[4] - u[3]), c("Colorado", "Sheep\nRiver", "Sheep\nMountain", "Dawson"), col = cols, xpd = NA)

legend(5, u[4], pch = c(1,19,8), pt.lwd= 1.5, col = 1, c("Common temp.\n common change", "Local temp.\ncommon change", "Local temp.\nlocal change"), xpd = NA, y.intersp = 1.5, bty = "n")

mtext(side = 3, outer = TRUE, c("Low emissions", "High emissions")[j - 1], line = -1, font = 2)

#---------------------------------------------
# Low and high together
quartz(width = 7, height = 3, pointsize = 10)
par(mfrow = c(1,2), mar = c(5,3,2,1), oma = c(0, 3, 0, 8))

for(j in 2:3){
	plot(1:4-0.1, AUC.all.percent[1, , j], ylim = extendrange(AUC.all.percent[, , j], f = 0.1), xlim = c(0.7, 4.3), col = cols, lwd = 2, xaxt = "n", xlab = "", las = 1, ylab = "")
	polygon(x = c(0,0,5,5), y = c(range(AUC.all.percent[3, , j]), rev(range(AUC.all.percent[3, , j]))), col = "#00000010", border =NA)
	points(1:4, AUC.all.percent[2, , j], col = cols, pch = 19)
	points(1:4 + 0.1, AUC.all.percent[3, , j], col = cols, pch = 8, lwd = 1.5)
	if(j == 2) mtext(side = 2, "Percent change in AUC\n under high emissions scenario", line = 3.5)
	
	u <- par('usr')
	abline(v = seq(1.5, 3.5, 1))
	axis(side = 1, at = c(1:4),, labels = FALSE)
	text(1:4, u[3] - 0.2*(u[4] - u[3]), c("Colorado", "Sheep\nRiver", "Sheep\nMountain", "Dawson"), col = cols, xpd = NA, srt = 90)
	mtext(side = 3, line = 0.5, paste0(letters[j-1], ") ", c("Low emissions", "High emissions")[j - 1]), adj = 0)
}

legend(4.7, u[4], pch = c(1,19,8), pt.lwd= 1.5, col = 1, c("Common temp.\n common change", "Local temp.\ncommon change", "Local temp.\nlocal change"), xpd = NA, y.intersp = 1.5, bty = "n")

################################################################################
# Look at transmission "window"
################################################################################

# Create array to store output
window.all <- array(NA, 
								 dim = c(3, 4, 3, 3), 
								 dimnames = list(c("scenario1", "scenario2", "scneario3"), LL[1:4], c("current", "rcp2.6", 'rcp8.5'), c("start", "end", "length")))

for(s in 1:3){ # for each scenario
	for(i in 1:4){ # For four locations, use parameters from those experiments
	for(j in 1:3){ # for different climate change scenarios
		window.all[s, i, j, 1:2] <- as.numeric(strftime(range(xDate[which(R0.all[s, , i, j] > 10^-10)]), format = "%j"))
		window.all[s, i, j, 3] <- diff(window.all[s, i, j, 1:2])
	}}}

# Plot length of window over different scenarios
quartz()
par(mfrow = c(3,1), mar = c(4,8,2,10))
for(s in 1:3){ # for each scenario
	plot(c(as.Date("2000-01-01"), as.Date("2000-12-31")), c(1,4), "n", xlab= "", ylab = "", yaxt = "n", ylim = c(0.5, 4.5), bty = "n")
	abline(v = c(as.Date(paste(2000, c(1:12), 1, sep = "-")), as.Date("2001-01-01")), lty = 3)
	axis(side = 2, at = c(1:4),c("Colorado", "Alberta", "South Yukon", "North Yukon")[1:4], las = 1)
	for(i in 1:4){
		segments(x0 = as.Date(paste(2000, window.all[s, i, 1, 1], sep = "-"), format = "%Y-%j"),
						 x1 = as.Date(paste(2000, window.all[s, i, 1, 2], sep = "-"), format = "%Y-%j"),
						 y0 = i, y1 = i, lwd = 16, col = grey(0.8))
		
		# Low & high emissions
		for(j in 2:3){
			segments(x0 = as.Date(paste(2000, window.all[s, i, j, 1], sep = "-"), format = "%Y-%j"),
						 x1 = as.Date(paste(2000, window.all[s, i, j, 2], sep = "-"), format = "%Y-%j"),
						 y0 = i + c(-0.1, 0.1)[j-1], y1 = i + c(-0.1, 0.1)[j-1], lwd = 3, col = c(4,2)[j-1], xpd = NA)
		} # end j
		
	} # end i
	mtext(side = 3, line = 0.5, adj = 0, paste0(letters[s], ") ", c("Common temp., common change", "Local temp., common change", "Local temp., local change")[s]))
	
	if(s == 1) legend(as.Date("2001-02-01"), 4, lwd = c(12, 3, 3), col = c(grey(0.8), 4, 2), c("Current", "Low emissions", "High emissions"), xpd = NA, bty = "n")
}
mtext(side = 1, expression(paste("Development window when ", D[L] > 0)), line = 3)
mtext(side = 2, outer= TRUE, "Source location for larvae", line = -1.2)

#--------------
# Change in length of window with climate projections
# Percent change in length of window

percWindow <- array(NA, dim = c(3, 4, 2), dimnames = list(c("scenario1", "scenario2", "scneario3"), LL[1:4], c("rcp2.6", 'rcp8.5')))
for(s in 1:3){
	for(i in 1:4){
		for(j in 1:2){
			percWindow[s, i, j] <- (window.all[s, i, j + 1, 3] - window.all[s, i, 1, 3])/(365-window.all[s, i, 1, 3])
		}
	}
}

quartz()
par(mfrow = c(3,1))
for(s in 1:3){
	plot(percWindow[s, , 1], 1:4, col = 4, lwd = 2, xlim = range(percWindow), ylim = c(0.5, 4.5), ylab = "", xlab = "")
points(percWindow[s, , 2], 1:4, col = 2, lwd = 2, cex = 1.5)	
}

#-------------------
# Plotting R0/C 
xDate <- as.Date(paste(2000, c(1:365), sep = "-"), format = "%Y-%j")

par(mfrow = c(3,1))
for(s in 1:3){
	plot(xDate, R0.all[s, , i, 1], "n", ylim = range(R0.all), las = 1, ylab = "")
	for(i in 1:4){
		polygon(x = c(xDate, rev(xDate)), y = c(R0.all[s, , i, 2], rev(R0.all[s, , i, 1])), col = paste0(cols[i], 30), border = NA)
		polygon(x = c(xDate, rev(xDate)), y = c(R0.all[s, , i, 3], rev(R0.all[s, , i, 1])), col = paste0(cols[i], 30), border = NA)
	}
}
