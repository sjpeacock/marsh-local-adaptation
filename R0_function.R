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

# Function to calculate R0C for different local climates and parasite parameters
# by simulating parasite dynamics under seasonally fluctuating temperature

calc_R0C <- function(
		climate, #  vector of daily temperatures, length 365
		parameters # named list of parameters mu and tau for temperaturese T.all = seq(-30, 40, 0.1)
		){
	
	DOY <- c(1:365)
	nYears <- 5
	
	# Join posterior draws are for T.all, so need to identify which
	# parameter to extract for each DOY based on tempLocChange
	tempInterval <- findInterval(climate, T.all)
	
	muL_instant <- as.numeric(parameters$mu[tempInterval])
	tauL_instant <- as.numeric(round(parameters$tau[tempInterval]))
	
	#-------------------
	# Numerically integrate tau
	#-------------------
	tau3 <- c(tauL_instant, tauL_instant, tauL_instant) # Pasted together for integration
	
	tauTilde <- numeric(length(DOY)) # How long has it taken an L3 at time t to develop from egg?
	tauHat <- numeric(length(DOY)) # How long does it take an egg shed at time t to develop to L3?
	for(j in DOY){
		
		# Tilde: How long has it taken an L3 at time t to develop from egg?
		dum <- as.numeric(cumsum(1/tau3[c((j + 365):1)]))
		tauTilde[j] <- which(dum > 1)[1]
		
		# Hat: How long does it take an egg shed at time t to develop to L3?
		boop <- as.numeric(cumsum(1/tau3[c((j + 365):length(tau3))]))
		tauHat[j] <- which(boop > 1)[1]
	}
	
	
	#-------------------
	# Extract instantaneous mortality and calulate DL
	#-------------------
	mu3 <- c(muL_instant, muL_instant, muL_instant)
	
	D_tilde <- numeric(length(DOY))
	D_hat <- numeric(length(DOY))
	for(j in DOY.nYears){
		
		# Tilde: What is the survival of an L3 at time t?
	D_tilde[j] <- exp(-sum(mu3[(j + 365 - tauTilde[j]):(j + 365)]))
		
		# What is the survival to L3 of an egg born at time t?
		D_hat[j] <- exp(-sum(mu3[(365 + j + 1):(365 + j + tauHat[j])]))
		
	}
	
	#-------------------
	# Calculate R0
	#-------------------
	R0 <- numeric(length(DOY))
	
	for(j in DOY){
		
		y <- array(0, dim = c(nYears * 365, 3), dimnames = list(NULL, c("L", "P", "R")))
		
		# Start with one egg at time t 
		# When is that egg a larvae? And what is the chance it survived?
		
		y[365 + j + tauHat[j], ] <- c(D_hat[j] * 1, 0, 0) # Initialize
		
		
		for(jj in (365 + j + tauHat[j]) : (365 * nYears - 1)){ # for the remaining time steps
			
			dy <- c(
				# Change in infectious larvae
				L = as.numeric(- muL_instant[DOY.nYears[jj]] * y[jj, "L"] - rho * H * y[jj, "L"] ),
				# Change in adult reproductive parasites in host
				P =  as.numeric(rho * H * y[jj - tauP, "L"] - (muP + bH) * y[jj, "P"] - alphaH * H * (y[jj, "P"]/H + (y[jj, "P"])^2/(H^2) * (kNB + 1)/kNB)),
				# Track reproductive output
				R = as.numeric(lambda  * y[jj, "P"]))#* D_tilde[DOY.nYears[i]] # Do not mutliply by DL
			
			y[jj + 1, ] <- pmax(0, y[jj, ] + dy)
			
		} # end jj
		
		R0[j] <- y[365 + j + 365*3, "R"]
		
			} # end day j
	
	R0C <- R0/(lambda * Dp/(alphaH + bH + muP))
	
	return(R0C)

} # end function