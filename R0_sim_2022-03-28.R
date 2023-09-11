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
# March 4, 2022
#
###############################################################################

###############################################################################
# Calculating R0(t) by simulating the set of parasite equations
# rather than R0(T) based on a constant temperature
# See Molnar et al. (2013 Ecology Letters)
#
# Across three scenarios, climate change, and locations
###############################################################################


library(dclone)
library(gplots)

# Colour for locations
cols <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')
params <- c("Colorado", "North Yukon")

# Tau and mu parameter estimates from -30 to 40 *C for four locations:
parPred <- readRDS("estimates/parPred.rds") 
# These estimates are based on join posterior draws and account for correlation
# among different hyperparameters (e.g., T_H and E_H)

# Temperature vector
T.all <- seq(-30, 40, 0.1)

# Temperature profiles for each location 
tempAll <- readRDS(file = "tempAll.rds")


# Set date variable
DOY <- c(1:365)
xDate <- as.Date(paste("2100", DOY, sep = "-"), format = "%Y-%j")

# All together
pdf(width = 3.2, height = 6, pointsize = 10, file = "figures/Final/tempAll.pdf")
par(mfrow = c(4,1), mar = c(3, 4.5, 2, 1), oma = c(0, 0, 3, 0))
for(i in 1:4){
	plot(xDate, tempAll[i, 1, 1, ], "l", ylim = range(tempAll[, , , ]), las = 1, ylab = "", xlab = "", xaxs = "i")
	abline(h = 0, col = grey(0.8))
	for(r in 1:3){
		for(k in 1:3){
			lines(xDate, tempAll[i, r, k, ], lty = k, col = c(1,4,2)[r], lwd = 0.8)
		}
	}
	mtext(side = 3, line = 0.5, adj = 0, paste0(letters[i], ") ", locNames[i]), cex = 0.8)
	
	if(i == 1) legend(xDate[90], 55, ncol = 2, col= c(1,4,2,1,1,1), lty = c(1,1,1,1,2,3), legend = c("historical", "RCP 2.6", "RCP 8.5", "Scenario 1", "Scenario 2", "Scenario 3"), xpd = NA, bty = "n", cex = 0.8)
}
mtext(side = 2, outer = TRUE, line = -1.5, expression(paste("Temperature (", degree, "C)")), cex = 0.8)

dev.off()

###############################################################################
# Simulate R0 over locations and scenarios
###############################################################################

#------------------------------------------------------------------------------
# Other parameters needed to simulate model
#------------------------------------------------------------------------------

kNB <- 1
lambda <- 1 # Shedding rate
rho <- 0.0001
H <- 0.01/rho # Assume rho*H = 0.01 as in paper
tauP <- 20 # take 20 days to develop from larvae to reproductive
muP <- 1/100 # Can live 100 days in host?
Dp <- exp(-muP * tauP)
alphaH <- 0.01 # Parasite-induced host mortality
bH <- 1/(7*365) # = 0.000391 Assume sheep live 7 years

#------------------------------------------------------------------------------
# Setup simulation
#------------------------------------------------------------------------------

nYears <- 5
DOY.nYears <- rep(DOY, nYears)

output <- list(
	tauL_instant = array(
		NA, 
		dim = c(2, 3, 3, 365), 
		dimnames = list(
			c("Colorado", "NorthYukon"), # Parameters to use 
			c("current", "rcp26", "rcp85"), 
			c("Scenario1", "Scenario2", "Scenario3"), # For North Yukon climate ONLY
			NULL)),
	muL_instant = array(NA, dim = c(2, 3, 3, 365), dimnames = list(c("Colorado", "NorthYukon"), c("current", "rcp26", "rcp85"), c("Scenario1", "Scenario2", "Scenario3"), NULL)),
	tauHat_integrated =  array(NA, dim = c(2, 3, 3, 365), dimnames = list(c("Colorado", "NorthYukon"), c("current", "rcp26", "rcp85"), c("Scenario1", "Scenario2", "Scenario3"), NULL)),
	tauTilde_integrated =  array(NA, dim = c(2, 3, 3, 365), dimnames = list(c("Colorado", "NorthYukon"), c("current", "rcp26", "rcp85"), c("Scenario1", "Scenario2", "Scenario3"), NULL)),
	DLHat_integrated =  array(NA, dim = c(2, 3, 3, 365), dimnames = list(c("Colorado", "NorthYukon"), c("current", "rcp26", "rcp85"), c("Scenario1", "Scenario2", "Scenario3"), NULL)),
	DLTilde_integrated =  array(NA, dim = c(2, 3, 3, 365), dimnames = list(c("Colorado", "NorthYukon"), c("current", "rcp26", "rcp85"), c("Scenario1", "Scenario2", "Scenario3"), NULL)),
	R0_integrated =  array(NA, dim = c(2, 3, 3, 365), dimnames = list(c("Colorado", "NorthYukon"), c("current", "rcp26", "rcp85"), c("Scenario1", "Scenario2", "Scenario3"), NULL)))


#------------------------------------------------------------------------------
# Simulate R0
#------------------------------------------------------------------------------

for(i in 1:2){ # For each location
	for(r in 1:3){ # For each climate change scenario
		for(k in c(1,2,3)){ # For each local scenario
			# Join posterior draws are for T.all, so need to identify which
			# parameter to extract for each DOY based on tempLocChange
			tempInterval <- findInterval(tempAll[4, r, k, ], T.all) # Use just North Yukon
			
			muL_instant <- as.numeric(parPred$mu[c(1,4)[i], tempInterval, 'mean'])
			output$muL_instant[i, r, k,  ] <- muL_instant
			
			tauL_instant <- as.numeric(round(parPred$tau[c(1,4)[i], tempInterval, 'mean']))
			output$tauL_instant[i, r, k, ] <- tauL_instant
			
			#-------------------
			# Numerically integrate tau
			#-------------------
			tau3 <- c(tauL_instant, tauL_instant, tauL_instant) # Pasted together for integration
			
			# tauTilde <- numeric(length(DOY)) # How long has it taken an L3 at time t to develop from egg?
			# tauHat <- numeric(length(DOY)) # How long does it take an egg shed at time t to develop to L3?
			for(j in DOY){
				
				# Tilde: How long has it taken an L3 at time t to develop from egg?
				dum <- as.numeric(cumsum(1/tau3[c((j + 365):1)]))
				output$tauTilde_integrated[i, r, k, j] <- which(dum > 1)[1]
				
				# Hat: How long does it take an egg shed at time t to develop to L3?
				boop <- as.numeric(cumsum(1/tau3[c((j + 365):length(tau3))]))
				output$tauHat_integrated[i, r, k, j] <- which(boop > 1)[1]
			}
			
			
			#-------------------
			# Extract instantaneous mortality and calulate DL
			#-------------------
			mu3 <- c(muL_instant, muL_instant, muL_instant)
			
			# D_tilde <- numeric(length(DOY))
			# D_hat <- numeric(length(DOY))
			for(j in DOY.nYears){
				
				# Tilde: What is the survival of an L3 at time t?
				output$DLTilde_integrated[i, r, k, j] <- exp(-sum(mu3[(j + 365 - output$tauTilde_integrated[i, r, k, j]):(j + 365)]))
				
				# What is the survival to L3 of an egg born at time t?
				output$DLHat_integrated[i, r, k, j] <- exp(-sum(mu3[(365 + j + 1):(365 + j + output$tauHat_integrated[i, r, k, j])]))
				
			}
			
			#-------------------
			# Calculate R0
			#-------------------
			for(j in DOY){
				
				y <- array(0, dim = c(nYears * 365, 3), dimnames = list(NULL, c("L", "P", "R")))
				
				# Start with one egg at time t 
				# When is that egg a larvae? And what is the chance it survived?
				
				y[365 + j + output$tauHat_integrated[i, r, k, j], ] <- c(output$DLHat_integrated[i, r, k, j] * 1, 0, 0) # Initialize
				
				
				for(jj in (365 + j + output$tauHat_integrated[i, r, k, j]) : (365 * nYears - 1)){ # for the remaining time steps
					
					dy <- c(
						# Change in infectious larvae
						L = as.numeric(- muL_instant[DOY.nYears[jj]] * y[jj, "L"] - rho * H * y[jj, "L"] ),
						# Change in adult reproductive parasites in host
						P =  as.numeric(rho * H * y[jj - tauP, "L"] - (muP + bH) * y[jj, "P"] - alphaH * H * (y[jj, "P"]/H + (y[jj, "P"])^2/(H^2) * (kNB + 1)/kNB)),
						# Track reproductive output
						R = as.numeric(lambda  * y[jj, "P"]))#* D_tilde[DOY.nYears[i]] # Do not mutliply by DL
					
					y[jj + 1, ] <- pmax(0, y[jj, ] + dy)
					
				} # end jj
				
				output$R0_integrated[i, r, k, j] <- y[365 + j + 365*3, "R"]
				
				
				# #-------------------
				# # Plot outcome
				# #-------------------
				# x <- c(DOY, 365 + DOY, 365*2 + DOY)
				# par(mfrow = c(3,1), mar = c(2,4,1,1), oma = c(2,0,2,0))
				# plot(x, y[365 + x, 'L'], "l", col = "#00B050", ylab = "L", lwd = 2, las = 1)
				# abline(v = j, lty = 2)
				# plot(x, y[365 + x, 'P'], "l", col = "#DF8344", ylab = "P", lwd = 2, las = 1)
				# abline(v = j, lty = 2)
				# plot(x, y[365 + x, 'R'], "l", col = "#38528B", ylab = "R0", lwd = 2, las = 1)
				# abline(v = j, lty = 2)
				
				
			} # end day j
			
			# plot(xDate, output$R0_integrated[i, r, k, ], "l", col = "#38528B", ylab = "R0", lwd = 2, las = 1, xaxs = "i")
			# axis(side = 3, at = as.Date(paste(seq(0, 365, 60), 2100, sep = "-"), format = "%j-%Y"), seq(0, 365, 60))
			
		} # end scenario k
	} # emissions r
} # end location i


# Tau tilde vs tau hat
plot(DOY, output$tauHat_integrated[1, 1, 1, ], "l", col = "purple")
lines(DOY, output$tauTilde_integrated[1, 1, 1, ], col = 2)
#------------------------------------------------------------------------------
# Calculate R0/C and RO_instant
#------------------------------------------------------------------------------

output$R0C <- output$R0_integrated /(lambda * Dp/(alphaH + bH + muP))

output$R0_instant <- exp(-output$muL_instant * output$tauL_instant) * rho * H /(output$muL_instant + rho * H)

#------------------------------------------------------------------------------
# Take cumulative sum of R0C (i.e., the AUC)
#------------------------------------------------------------------------------

output$cumR0C <- array(NA, dim = c(2, 3, 3), dimnames = list(c("Colorado", "NorthYukon"), c("current", "rcp26", "rcp85"), c("Scenario1", "Scenario2", "Scenario3")))

for(i in 1:2){ # For each location
	for(r in 1:3){ # For each climate change scenario
		for(k in 1:3){ # For each local scenario
			
			output$cumR0C[i, r, k] <- sum(output$R0C[i, r, k, ])
			
		} # end scenario k
	} # emissions r
} # end location i

###############################################################################
# Filling in AUC table
###############################################################################

plot(xDate, output$R0C[1, 1, 1, ], ylim = range(output$R0C), "l")
lines(xDate, output$R0C[1, 2, 1, ], col = 4)
lines(xDate, output$R0C[1, 3, 1, ], col = 2)

i <- 2
k <- 1
# quartz(width = 3, height = 2.5, pointsize = 10)
par(mar = c(2,2,1,1))
plot(xDate, output$R0C[i, 3, k, ]-output$R0C[i, 1, k, ], ylab = "Change in R0C", "l", col = 2, ylim = range(output$R0C[, 3, , ]-output$R0C[, 1, , ]))
lines(xDate, output$R0C[i, 2, k, ]-output$R0C[i, 1, k, ], col = 4)
abline(h = 0, lty = 2)

plot(xDate, output$R0C[1, 3, 1, ]-output$R0C[1, 1, 1, ], ylab = "Change in R0C", "l", col = 2, ylim = range(output$R0C[, 3, , ]-output$R0C[, 1, , ]), lwd = 2, lty = 2)
lines(xDate, output$R0C[2, 3, 3, ]-output$R0C[2, 1, 3, ], lwd = 2, col = 2)

lines(xDate, output$R0C[2, 3, 1, ]-output$R0C[2, 1, 1, ], lty = 3, col = 2)

lines(xDate, output$R0C[1, 3, 3, ]-output$R0C[1, 1, 3, ], lty = 3, col = 2)

(sum(output$R0C[2, 3, 3, ])-sum(output$R0C[2, 1, 3, ]))/sum(output$R0C[2, 1, 3, ])

(sum(output$R0C[2, 3, 3, ])-sum(output$R0C[2, 1, 3, ]))/sum(output$R0C[2, 1, 3, ])


###############################################################################
###############################################################################

# Compare Colorado climate and North Yukon parameters

# Compare climate (left) using local parameters
#   Colorado vs Colorado
#   Colorado vs Alberta
#   Colorado vs South Yukon
#   Colorado vs North Yukon

# Compare parameters using local climate
#   Colorado vs Colorado
#   Colorado vs Alberta
#   Colorado vs South Yukon
#   Colorado vs North Yukon


#------------------------------------------------------------------------------
# Source function to integrate and calculate R0
source("R0_function.R")

# Setup array
R0C <- list(
	CompareClimate = array(data = NA, dim = c(2, 4, 365)),
	CompareParameters = array(data = NA, dim = c(2, 4, 365))
)

# tempAll[location, rcp, scenario1-3, DOY]
# r <- 1
# s <- 3

# Compare climate using local parameter
for(i in 1:4){
	R0C[[1]][1, i, ] <- calc_R0C(
		climate = tempAll[1, 1, 3, ], # Colorado climate
		parameters = list(# Local parameters
			mu = as.numeric(parPred$mu[i, , 'mean']),
			tau = as.numeric(round(parPred$tau[i, , 'mean']))
			)) 

	R0C[[1]][2, i, ] <- calc_R0C(
		climate = tempAll[i, 1, 3, ], # Local climate
		parameters = list(# Local parameters
			mu = as.numeric(parPred$mu[i, , 'mean']),
			tau = as.numeric(round(parPred$tau[i, , 'mean']))
		)) 
}

# Compare parameters using local climate
for(i in 1:4){
	R0C[[2]][1, i, ] <- calc_R0C(
		climate = tempAll[i, 1, 3, ], # Local climate
		parameters = list(# Colorado parameters
			mu = as.numeric(parPred$mu[1, , 'mean']),
			tau = as.numeric(round(parPred$tau[1, , 'mean']))
		)) 
	
	R0C[[2]][2, i, ] <- calc_R0C(
		climate = tempAll[i, 1, 3, ], # Local climate
		parameters = list(# Local parameters
			mu = as.numeric(parPred$mu[i, , 'mean']),
			tau = as.numeric(round(parPred$tau[i, , 'mean']))
		)) 
}

#------------------------------------------------------------------------------
# Comparing Scenario 1 & 3 under current conditions
#------------------------------------------------------------------------------
xTicks <- as.Date(paste("01", c(1:12), "2100", sep = "-"), format = "%d-%m-%Y")

quartz(width = 6.3, height = 7, pointsize = 10)
par(mfcol = c(4,2), mar = c(0,0,0,0), oma = c(5, 6, 3, 1))

for(j in 1:2){
for(i in 1:4){
	plot(xDate, R0C[[j]][1, i, ], "n", ylim = c(0.3, 1.2), xlab = "", xaxt = "n", xaxs = "i", yaxt = "n")
	if(j == 1) axis(side = 2, las = 1) else axis(side = 2, labels = FALSE)
	abline(h = seq(0.3, 1.2, 0.1), lwd = 0.5, col = grey(0.8))
	abline(v = as.Date("2100-09-18"), col = 2, lwd = 2)
	lines(xDate, R0C[[j]][1, i, ], lty = 2)
	lines(xDate, R0C[[j]][2, i, ], lwd = 1.5)
	axis(side = 1, at = xTicks[seq(1, 12, 2)], labels = FALSE)
	axis(side = 1, at = xTicks, labels = FALSE, tck = -0.02)
	mtext(side = 3, line = -2.5, paste0("  ", letters[(j-1)*4 + i], ") ", locNames[i]), adj=0)
	
	if(i==1 & j == 1) mtext(side = 3, "Local parameters", line = 1.5)
	if(i==1 & j == 2) mtext(side = 3, "Local climate", line = 1.5)
	
}

axis(side = 1, at = xTicks[seq(1, 12, 2)], labels = month.abb[seq(1, 12, 2)])

if(j == 1){
	legend("bottomleft", lty = c(2,1), lwd = c(1,1.5), c("Colorado", "Local"), title = "Climate", bg = "white")
mtext(side = 2, line = 3.5, outer = TRUE, "Relative parasite performance (R0/C)")
} else {
	legend("bottomleft", lty = c(2,1), lwd = c(1,1.5), c("Colorado", "Local"), title = "Parameters", bg = "white")
}
} # end j


###############################################################################
###############################################################################
# April 11, 2023
# Compare all parameters and all climates under all climate change scnearios
# (current, rcp45, rcp85)
###############################################################################
###############################################################################

locNames <- c("Colorado", "Alberta", "South Yukon", "North Yukon")

#------------------------------------------------------------------------------
# Source function to integrate and calculate R0
source("R0_function.R")

# Setup array
R0C <- array(data = NA,
						 dim = c(4, 4, 3, 365),
						 dimnames = list(
						 	c(paste0("Climate:", locNames)),
						 	c(paste0("Parameters:", locNames)),
						  c("current", "rcp45", "rcp85"),
						 	NULL))


for(i in 1:4){ # For each climate
	for(j in 1:4){ # For each parameters
		for(r in 1:3){ # For each climate change scenario
			R0C[i, j, r, ] <- calc_R0C(
				climate = tempAll[i, r, 3, ], # Colorado climate
				parameters = list(# Local parameters
					mu = as.numeric(parPred$mu[j, , 'mean']),
					tau = as.numeric(round(parPred$tau[j, , 'mean']))
				)) 
}}}


# Write to csv
R0C_out <- data.frame(
	Climate = rep(locNames, each = 4*3*365),
	Parameters = rep(rep(locNames, each = 3*365), 4),
	ClimateChange = rep(rep(c("current", "rcp45", "rcp85"), each = 365), 4*4),
	DOY = rep(c(1:365), 4*4*3),
	R0C = rep(NA, 4*4*3*365))
for(i in 1:4){ # For each climate
	for(j in 1:4){ # For each parameters
		for(r in 1:3){ # For each climate change scenario
			R0C_out[which(R0C_out$Climate == locNames[i] &R0C_out$Parameters == locNames[j] & R0C_out$ClimateChange == c("current", "rcp45", "rcp85")[r]), "R0C"] <- R0C[i, j, r, ]
		}}}	

write.csv(R0C_out, file = "output/R0C_all_combos_2023-04-11.csv")

# Calculate AUC
AUC_out <- data.frame(
	Climate = rep(locNames, each = 4*3),
	Parameters = rep(rep(locNames, each = 3), 4),
	ClimateChange = rep(c("current", "rcp45", "rcp85"), 4*4),
	AUC = rep(NA, 4*4*3))
for(i in 1:4){ # For each climate
	for(j in 1:4){ # For each parameters
		for(r in 1:3){ # For each climate change scenario
			AUC_out[which(AUC_out$Climate == locNames[i] &AUC_out$Parameters == locNames[j] & AUC_out$ClimateChange == c("current", "rcp45", "rcp85")[r]), "AUC"] <- sum(R0C[i, j, r, ])
		}}}	

write.csv(AUC_out, file = "output/AUC_all_combos_2023-04-11.csv", row.names = FALSE)

# Calculate % Change
perc <- data.frame(
	Climate = rep(rep(locNames, each = 4), 2),
	Parameters = rep(locNames, 4*2),
	ClimateChange = rep(c("rcp45", "rcp85"), each = 4*4),
	percChange = rep(NA, 4*4*2)
)
for(i in 1:4){ # For each climate
	for(j in 1:4){ # For each parameters
		for(r in 1:2){ # For each climate change scenario
			perc[which(perc$Climate == locNames[i] & perc$Parameters == locNames[j] & perc$ClimateChange == c("rcp45", "rcp85")[r]), "percChange"] <- 
				(AUC_out$AUC[which(AUC_out$Climate == locNames[i] &AUC_out$Parameters == locNames[j] & AUC_out$ClimateChange == c("rcp45", "rcp85")[r])] - AUC_out$AUC[which(AUC_out$Climate == locNames[i] & AUC_out$Parameters == locNames[j] & AUC_out$ClimateChange == "current")])/AUC_out$AUC[which(AUC_out$Climate == locNames[i] & AUC_out$Parameters == locNames[j] & AUC_out$ClimateChange == "current")] * 100
		}}}	
	
write.csv(perc, file = "output/percChange_all_combos_2023-04-11.csv", row.names = FALSE)
