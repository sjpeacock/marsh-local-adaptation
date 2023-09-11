# Thermal tolerance and local adaptation to freezing temperatures in a key nematode of Arctic ungulates: 
# the case of Marshallagia marshalli 

# File with model functions for fitting in JAGS or R optim
# Note: fits with optim() have not been tested.

# Author: Stephanie Peacock <stephanie.j.peacock@gmail.com>
# Data owner: O. Alejandro Aleuy <oaleuy@ucalgary.ca>


###############################################################################
# Metabolic functions
###############################################################################
predictMTE <- function(Temp, params, fn="A"){
	if(length(which(names(params) == 'mu0')) == 0){
		a <- params['tau0']
		Efront <- 1
	} else {
		a <- params['mu0']
		Efront <- -1
	}
	if(fn=="A"){
		parEst <- a * exp(Efront * params['E']/(8.62*10^-5) * (1/(Temp + 273.15) - 1/(15 + 273.15)))
	
		} else if (fn=="SSU"){
			parEst <- a * exp(Efront * params['E']/(8.62*10^-5) * (1/(Temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['Eh']/(8.62*10^-5) * (-1/(Temp + 273.15) + 1/(params['Th'] + 273.15))))
	
		} else if (fn=="SSUL"){
			parEst <- a * exp(Efront * params['E']/(8.62*10^-5) * (1/(Temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['El']/(8.62*10^-5) * (1/(Temp + 273.15) - 1/(params['Tl'] + 273.15))) + exp(params['Eh']/(8.62*10^-5) * (-1/(Temp + 273.15) + 1/(params['Th'] + 273.15))))
	
		} else if (fn=="SSL"){
			parEst <- a * exp(Efront * params['E']/(8.62*10^-5) * (1/(Temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['El']/(8.62*10^-5) * (1/(Temp + 273.15) - 1/(params['Tl'] + 273.15))))
		
			} else {
		stop("Function not defined.")
	}
	return(parEst)
}

###############################################################################
# JAGS model for survival
###############################################################################

modelSurv <- function(){
	
	# Priors on parameters
	# sigma_v ~ dunif(0,100)
	for(i in 1:nT){
		for(j in 1:nL){
			mu[i,j] ~ dlnorm(log(0.1), 1) # First guess
			}}
	
	# Model
	# # Lack of fit of incubator k 
	# for(j in 1:nj) {
	# 	v[j] ~ dnorm(0, sigma_v^(-2))
	# }

	# Calculate probability for each data point
	for(z in 1:n) {
		# p[z] <- (1 - censored[z])*(exp(- exp(v[incubator[z]]) * mu[temp[z], loc[z]] * time_back[z]) - exp(- exp(v[incubator[z]]) * mu[temp[z], loc[z]] * time[z])) + censored[z]*(exp(- exp(v[incubator[z]]) * mu[temp[z], loc[z]]*time_back[z]))
		p[z] <- (1 - censored[z])*(exp(- mu[temp[z], loc[z]] * time_back[z]) - exp(- mu[temp[z], loc[z]] * time[z])) + censored[z]*(exp(- mu[temp[z], loc[z]]*time_back[z]))
		}
	
	# Likleihood
	for(z in 1:n) {
		y[z] ~ dbern(max(10^-10, p[z]))
	}


} # end model
	
	# Parameters to estimate: 
# tau_t, sigma_v

# censored
# parameters: sigma_v, sigma_e, tau (length nT)
# data: 
#		nT: number of temperature treatments
# 	nj: number of incubrators for each temperature
#		time: timestep at which development was observed (matrix nT, nj)
#		time_back: last timestep observed before development (matrix nT, nj)
# 	censored: was the data censored? 0 = no, 1 = yes (matrix nT, nj) 

###############################################################################
# JAGS model for development time
###############################################################################

modelSurvMTE.A <- function(){
	
	# Priors on parameters
	for(l in 1:nL){
		mu0[l] ~ dlnorm(log(30), 1)
		E[l] ~ dlnorm(log(0.65), 3) #plot(seq(0, 1.2, 0.01), dlnorm(seq(0, 1.2, 0.01), log(0.65), 3^(-2)), "l")
		}
	
	# MTE model for development time: eq (3a) from Molnar et al. (2013) 
	for(i in 1:nT){
		for(l in 1:nL){
			mu[i,l] <- mu0[l] * exp(-E[l]/(8.62*10^-5) * (1/(T.obs[i] + 273.15) - 1/(15 + 273.15)))
		}}
	
	# Calculate probability for each data point
	for(z in 1:n) {
		p[z] <- (1 - censored[z])*(exp(- mu[temp[z], loc[z]] * time_back[z]) - exp(- mu[temp[z], loc[z]] * time[z])) + censored[z]*(exp(- mu[temp[z], loc[z]]*time_back[z]))
	}
	
	# Likleihood
	for(z in 1:n) {
		y[z] ~ dbern(max(10^-10, p[z]))
	}
	
} # end MTE model

modelSurvMTE.SSU <- function(){
	
	# Priors on parameters
	for(l in 1:nL){
		mu0[l] ~ dlnorm(log(30), 1)
		E[l] ~ dlnorm(log(0.65), 3) #plot(seq(0, 1.2, 0.01), dlnorm(seq(0, 1.2, 0.01), log(0.65), 3^(-2)), "l")
		Eh[l] ~ dlnorm(log(3.65), 3)
		Th[l] ~ dnorm(25, 5^(-2))
		}
	
	# MTE model for development time: eq (3a) from Molnar et al. (2013) 
	for(i in 1:nT){
		for(l in 1:nL){
			mu[i,l] <- mu0[l] * exp(- E[l]/(8.62*10^-5) * (1/(T.obs[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(Eh[l]/(8.62*10^-5) * (-1/(T.obs[i] + 273.15) + 1/(Th[l] + 273.15))))
		}}
	
	# Model
	# Calculate probability for each data point
	for(z in 1:n) {
		p[z] <- (1 - censored[z])*(exp(- mu[temp[z], loc[z]] * time_back[z]) - exp(- mu[temp[z], loc[z]] * time[z])) + censored[z]*(exp(- mu[temp[z], loc[z]]*time_back[z]))
	}
	
	# Likleihood
	for(z in 1:n) {
		y[z] ~ dbern(max(10^-10, p[z]))
	}
	
} # end MTE model

# modelSurvMTE.SSL <- function(){
# 	
# 	# Priors on parameters
# 	for(l in 1:nL){
# 		mu0[l] ~ dlnorm(log(30), 1)
# 		E[l] ~ dlnorm(log(0.65), 3) #plot(seq(0, 1.2, 0.01), dlnorm(seq(0, 1.2, 0.01), log(0.65), 3^(-2)), "l")
# 		El[l] ~ dlnorm(log(3.65), 1)
# 		Tl[l] ~ dnorm(5, 5^(-2))
# 		}
# 	
# 	# MTE model for development time: eq (3a) from Molnar et al. (2013) 
# 	for(i in 1:nT){
# 		for(l in 1:nL){
# 			mu[i,l] <- mu0[l] * exp( - E[l]/(8.62*10^-5) * (1/(T.obs[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(El[l]/(8.62*10^-5) * (1/(T.obs[i] + 273.15) - 1/(Tl[l] + 273.15))))
# 		}}
# 	
# 	# Model
# 	# Calculate probability for each data point
# 	for(z in 1:n) {
# 		p[z] <- (1 - censored[z])*(exp(- mu[temp[z], loc[z]] * time_back[z]) - exp(- mu[temp[z], loc[z]] * time[z])) + censored[z]*(exp(- mu[temp[z], loc[z]]*time_back[z]))
# 	}
# 	
# 	# Likleihood
# 	for(z in 1:n) {
# 		y[z] ~ dbern(max(10^-10, p[z]))
# 	}
# 	
# } # end MTE model

modelSurvMTE.SSUL <- function(){
	
	# Priors on parameters
	for(j in 1:nL){
		mu0[j] ~ dlnorm(log(30), 1)
		E[j] ~ dlnorm(log(0.65), 3) #plot(seq(0, 100, 0.01), dlnorm(seq(0, 100, 0.01), log(30), 1^(-2)), "l")
		El[j] ~ dlnorm(log(3.65), 3)
		Eh[j] ~ dlnorm(log(3.65), 3)
		Tl[j] ~ dnorm(5, 5^(-2))
		Th[j] ~ dnorm(25, 5^(-2))
	}
	
	# MTE model for development time: eq (3a) from Molnar et al. (2013) 
	for(i in 1:nT){
		for(j in 1:nL){
			mu[i,j] <- mu0[j] * exp( - E[j]/(8.62*10^-5) * (1/(T.obs[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(El[j]/(8.62*10^-5) * (1/(T.obs[i] + 273.15) - 1/(Tl[j] + 273.15))) + exp(Eh[j]/(8.62*10^-5) * (-1/(T.obs[i] + 273.15) + 1/(Th[j] + 273.15)))) 
		}}
	
	# Model
	# Calculate probability for each data point
	for(z in 1:n) {
		p[z] <- (1 - censored[z])*(exp(- mu[temp[z], loc[z]] * time_back[z]) - exp(- mu[temp[z], loc[z]] * time[z])) + censored[z]*(exp(- mu[temp[z], loc[z]]*time_back[z]))
	}
	# Likleihood
	for(z in 1:n) {
		y[z] ~ dbern(max(10^-10, p[z]))
	}
	
} # end MTE model
# ###############################################################################
# # logLikelihood version for optim() in R
# ###############################################################################
logLikSurv <- function(params, dat.list, mod){
	
	o <- names(params)
	
		
	mu <- matrix(NA, nrow = dat.list$nT, ncol = dat.list$nL)
	
	#----------------------------------
	if(mod == "IT"){
		
		for(i in 1:dat.list$nT){
			for(j in 1:dat.list$nL){
				mu[i,j] <- as.numeric(params[paste("mu[", i, ",", j, "]", sep="")])
			}}
	
	#----------------------------------
	} else if(mod == "SSU"){
	
		for(i in 1:dat.list$nT){
			for(j in 1:dat.list$nL){
				mu[i,j] <- as.numeric(params[paste("mu0[", j, "]", sep="")] * exp(- params[paste("E[", j, "]", sep="")]/(8.62*10^-5) * (1/(dat.list$T.obs[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(params[paste("Eh[", j, "]", sep="")]/(8.62*10^-5) * (-1/(dat.list$T.obs[i] + 273.15) + 1/(params[paste("Th[", j, "]", sep="")] + 273.15)))))
				
			}}
		
	#----------------------------------
	} else if(mod == "SSUL"){
		for(i in 1:dat.list$nT){
			for(j in 1:dat.list$nL){
				mu[i,j] <- as.numeric(params[paste("mu0[", j, "]", sep="")] * exp(- params[paste("E[", j, "]", sep="")]/(8.62*10^-5) * (1/(dat.list$T.obs[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(params[paste("El[", j, "]", sep="")]/(8.62*10^-5) * (1/(dat.list$T.obs[i] + 273.15) - 1/(params[paste("Tl[", j, "]", sep="")] + 273.15))) + exp(params[paste("Eh[", j, "]", sep="")]/(8.62*10^-5) * (-1/(dat.list$T.obs[i] + 273.15) + 1/(params[paste("Th[", j, "]", sep="")] + 273.15)))))
				
			}}
		#----------------------------------
	} else if(mod == "SSU_1L"){
		
		for(i in 1:dat.list$nT){
			for(j in 1:dat.list$nL){
				mu[i,j] <- as.numeric(params["mu0"] * exp(- params["E"]/(8.62*10^-5) * (1/(dat.list$T.obs[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(params["Eh"]/(8.62*10^-5) * (-1/(dat.list$T.obs[i] + 273.15) + 1/(params["Th"] + 273.15)))))
			}}
		#----------------------------------
	} else if(mod == "A"){
		
		for(i in 1:dat.list$nT){
			for(j in 1:dat.list$nL){
				mu[i,j] <- as.numeric(params[paste("mu0[", j, "]", sep="")] * exp(- params[paste("E[", j, "]", sep="")]/(8.62*10^-5) * (1/(dat.list$T.obs[i] + 273.15) - 1/(15 + 273.15))))
			}}
		
		#----------------------------------
	}
	
	# Model
	# Calculate probability for each data point
	p <- numeric(dat.list$n)
	for(z in 1:dat.list$n) {
		p[z] <- (1 - dat.list$censored[z])*(exp(- mu[dat.list$temp[z], dat.list$loc[z]] * dat.list$time_back[z]) - exp(- mu[dat.list$temp[z], dat.list$loc[z]] * dat.list$time[z])) + dat.list$censored[z]*(exp(- mu[dat.list$temp[z], dat.list$loc[z]]*dat.list$time_back[z]))
	}
	
	# Likleihood
	LLd <- numeric(dat.list$n)
	for(z in 1:dat.list$n) {
		LLd[z] <- dbinom(1, size = 1, prob = max(p[z], 10^-10), log=TRUE)
	}

	return(list(NLL = -sum(LLd), -LLd))
		}