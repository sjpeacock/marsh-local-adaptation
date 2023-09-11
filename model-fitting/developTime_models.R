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


# File with model functions for fitting in JAGS or R optim
# Note: fits with optim() have not been tested.


###############################################################################
# Metabolic functions
###############################################################################
predictMTE <- function(Temp, params, fn="A"){
	
	if(fn=="A"){
		tau <- params['tau0'] * exp(params['E']/(8.62*10^-5) * (1/(Temp + 273.15) - 1/(15 + 273.15)))
	
		} else if (fn=="SSU"){
		tau <- params['tau0'] * exp(params['E']/(8.62*10^-5) * (1/(Temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['Eh']/(8.62*10^-5) * (-1/(Temp + 273.15) + 1/(params['Th'] + 273.15))))
	
		} else if (fn=="SSUL"){
		tau <- params['tau0'] * exp(params['E']/(8.62*10^-5) * (1/(Temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['El']/(8.62*10^-5) * (1/(Temp + 273.15) - 1/(params['Tl'] + 273.15))) + exp(params['Eh']/(8.62*10^-5) * (-1/(Temp + 273.15) + 1/(params['Th'] + 273.15))))
	
		} else if (fn=="SSL"){
			tau <- params['tau0'] * exp(params['E']/(8.62*10^-5) * (1/(Temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['El']/(8.62*10^-5) * (1/(Temp + 273.15) - 1/(params['Tl'] + 273.15))))
		
			} else {
		stop("Function not defined.")
	}
	return(tau)
}

###############################################################################
# JAGS model for development time
###############################################################################

modelDevelop <- function(){
	
	# Priors on parameters
	sigma_v ~ dunif(0,100)
	sigma_e ~ dunif(0, 100)
	for(i in 1:nT){
		for(j in 1:nL){
			tau[i,j] ~ dexp(1/100) #Changed from 1/15 which gave wierd results for Sheep Mountain (few uncensored data)
			}}
	
	# Model
	# Lack of fit of incubator k (upsilon in Regniere et al. (2012))
	for(j in 1:nj) {
		v[j] ~ dnorm(1, sigma_v^(-2))
	}

	for(z in 1:n) {
		tpred[z] <- tau[temp[z], loc[z]]*v[incubator[z]] 
		# S[z] <- gamma[temp[z], loc[z]]^(tau[temp[z], loc[z]]*v[incubator[z]]) #Survival to L3
	}
	
	
	# Calculate epsilon for each data point, and probability.
	for(z in 1:n) {
		eps[z] <- log(time[z]/tpred[z]) #epsilon: error for that point
		eps_back[z] <- log(time_back[z]/tpred[z]) #epsilon minus 1
		p[z] <- (1 - censored[z])*(phi((eps[z] + 0.5*sigma_e^2)/sigma_e) - phi((eps_back[z] + 0.5*sigma_e^2)/sigma_e)) + censored[z]*(1 - phi((eps[z] + 0.5*sigma_e^2)/sigma_e))
		}
	
	# Likleihood
	
	for(z in 1:n) {
		y[z] ~ dbern(max(10^-10, p[z]))
	}


} # end model
	
	# Parameters to estimate: 
# tau_t, sigma_v, sigma_e

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

modelDevelopMTE.A <- function(){
	
	# Priors on parameters
	sigma_v ~ dunif(0,100)
	sigma_e ~ dunif(0, 100)
	for(l in 1:nL){
		tau0[l] ~ dlnorm(log(30), 1)
		E[l] ~ dlnorm(log(0.65), 3) #plot(seq(0, 1.2, 0.01), dlnorm(seq(0, 1.2, 0.01), log(0.65), 3^(-2)), "l")
		}
	
	# MTE model for development time: eq (3a) from Molnar et al. (2013) 
	for(i in 1:nT){
		for(l in 1:nL){
			tau[i,l] <- tau0[l] * exp(E[l]/(8.62*10^-5) * (1/(T.obs[i] + 273.15) - 1/(15 + 273.15)))
		}}
	
	# Model
	# Lack of fit of incubator k (upsilon in Regniere et al. (2012))
	for(j in 1:nj) {
		v[j] ~ dnorm(1, sigma_v^(-2))
		}
	
	for(z in 1:n) {
		tpred[z] <- tau[temp[z], loc[z]]*v[incubator[z]] 
		# S[z] <- gamma[temp[z], loc[z]]^(tau[temp[z], loc[z]]*v[incubator[z]]) #Survival to L3
	}
	
	
	# Calculate epsilon for each data point, and probability.
	for(z in 1:n) {
		eps[z] <- log(time[z]/tpred[z]) #epsilon: error for that point
		eps_back[z] <- log(time_back[z]/tpred[z]) #epsilon minus 1
		p[z] <- (1 - censored[z])*(phi((eps[z] + 0.5*sigma_e^2)/sigma_e) - phi((eps_back[z] + 0.5*sigma_e^2)/sigma_e)) + censored[z]*(1 - phi((eps[z] + 0.5*sigma_e^2)/sigma_e))
	}
	
	# Likleihood
	
	for(z in 1:n) {
		y[z] ~ dbern(max(10^-10, p[z]))
	}
	
} # end MTE model

modelDevelopMTE.SSU <- function(){
	
	# Priors on parameters
	sigma_v ~ dunif(0,100)
	sigma_e ~ dunif(0, 100)
	for(l in 1:nL){
		tau0[l] ~ dlnorm(log(30), 1)
		E[l] ~ dlnorm(log(0.65), 3) #plot(seq(0, 1.2, 0.01), dlnorm(seq(0, 1.2, 0.01), log(0.65), 3^(-2)), "l")
		Eh[l] ~ dlnorm(log(3.65), 3)
		Th[l] ~ dnorm(25, 5^(-2))
		}
	
	# MTE model for development time: eq (3a) from Molnar et al. (2013) 
	for(i in 1:nT){
		for(l in 1:nL){
			tau[i,l] <- tau0[l] * exp(E[l]/(8.62*10^-5) * (1/(T.obs[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(Eh[l]/(8.62*10^-5) * (-1/(T.obs[i] + 273.15) + 1/(Th[l] + 273.15)))) #exp(El[j]/(8.62*10^-5) *lll (1/(T.obs[i] + 273.15) - 1/(Tl[j] + 273.15))) +
		}}
	
	# Model
	# Lack of fit of incubator k (upsilon in Regniere et al. (2012))
	for(j in 1:nj) {
		v[j] ~ dnorm(1, sigma_v^(-2))
	}

	for(z in 1:n) {
		tpred[z] <- tau[temp[z], loc[z]]*v[incubator[z]] 
		# S[z] <- gamma[temp[z], loc[z]]^(tau[temp[z], loc[z]]*v[incubator[z]]) #Survival to L3
	}
	
	
	# Calculate epsilon for each data point, and probability.
	for(z in 1:n) {
		eps[z] <- log(time[z]/tpred[z]) #epsilon: error for that point
		eps_back[z] <- log(time_back[z]/tpred[z]) #epsilon minus 1
		p[z] <- (1 - censored[z])*(phi((eps[z] + 0.5*sigma_e^2)/sigma_e) - phi((eps_back[z] + 0.5*sigma_e^2)/sigma_e)) + censored[z]*(1 - phi((eps[z] + 0.5*sigma_e^2)/sigma_e))
	}
	
	# Likleihood
	
	for(z in 1:n) {
		y[z] ~ dbern(max(10^-10, p[z]))
	}
	
} # end MTE model

modelDevelopMTE.SSL <- function(){
	
	# Priors on parameters
	sigma_v ~ dunif(0, 100)
	sigma_e ~ dunif(0, 100)
	for(l in 1:nL){
		tau0[l] ~ dlnorm(log(30), 1)
		E[l] ~ dlnorm(log(0.65), 3) #plot(seq(0, 1.2, 0.01), dlnorm(seq(0, 1.2, 0.01), log(0.65), 3^(-2)), "l")
		El[l] ~ dlnorm(log(3.65), 1)
		Tl[l] ~ dnorm(5, 5^(-2))
		}
	
	# MTE model for development time: eq (3a) from Molnar et al. (2013) 
	for(i in 1:nT){
		for(l in 1:nL){
			tau[i,l] <- tau0[l] * exp(E[l]/(8.62*10^-5) * (1/(T.obs[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(El[l]/(8.62*10^-5) * (1/(T.obs[i] + 273.15) - 1/(Tl[l] + 273.15)))) # + exp(Eh[j]/(8.62*10^-5) * (-1/(T.obs[i] + 273.15) + 1/(Th[j] + 273.15)))
		}}
	
	# Model
	# Lack of fit of incubator k (upsilon in Regniere et al. (2012))
	for(j in 1:nj) {
		v[j] ~ dnorm(1, sigma_v^(-2))
	}
	
	for(z in 1:n) {
		tpred[z] <- tau[temp[z], loc[z]]*v[incubator[z]] 
		# S[z] <- gamma[temp[z], loc[z]]^(tau[temp[z], loc[z]]*v[incubator[z]]) #Survival to L3
	}
	
	
	# Calculate epsilon for each data point, and probability.
	for(z in 1:n) {
		eps[z] <- log(time[z]/tpred[z]) #epsilon: error for that point
		eps_back[z] <- log(time_back[z]/tpred[z]) #epsilon minus 1
		p[z] <- (1 - censored[z])*(phi((eps[z] + 0.5*sigma_e^2)/sigma_e) - phi((eps_back[z] + 0.5*sigma_e^2)/sigma_e)) + censored[z]*(1 - phi((eps[z] + 0.5*sigma_e^2)/sigma_e))
	}
	
	# Likleihood
	
	for(z in 1:n) {
		y[z] ~ dbern(max(10^-10, p[z]))
	}
	
} # end MTE model

modelDevelopMTE.SSUL <- function(){
	
	# Priors on parameters
	sigma_v ~ dunif(0,100)
	sigma_e ~ dunif(0,100)
	for(j in 1:nL){
		tau0[j] ~ dlnorm(log(30), 1)
		E[j] ~ dlnorm(log(0.65), 3) #plot(seq(0, 100, 0.01), dlnorm(seq(0, 100, 0.01), log(30), 1^(-2)), "l")
		El[j] ~ dlnorm(log(3.65), 3)
		Eh[j] ~ dlnorm(log(3.65), 3)
		Tl[j] ~ dnorm(5, 5^(-2))
		Th[j] ~ dnorm(25, 5^(-2))
	}
	
	# MTE model for development time: eq (3a) from Molnar et al. (2013) 
	for(i in 1:nT){
		for(j in 1:nL){
			tau[i,j] <- tau0[j] * exp(E[j]/(8.62*10^-5) * (1/(T.obs[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(El[j]/(8.62*10^-5) * (1/(T.obs[i] + 273.15) - 1/(Tl[j] + 273.15))) + exp(Eh[j]/(8.62*10^-5) * (-1/(T.obs[i] + 273.15) + 1/(Th[j] + 273.15)))) 
		}}
	
	# Model
	# Lack of fit of incubator k (upsilon in Regniere et al. (2012))
	for(j in 1:nj) {
		v[j] ~ dnorm(1, sigma_v^(-2))
	}
	
	for(z in 1:n) {
		tpred[z] <- tau[temp[z], loc[z]]*v[incubator[z]] 
		# S[z] <- gamma[temp[z], loc[z]]^(tau[temp[z], loc[z]]*v[incubator[z]]) #Survival to L3
	}
	
	# Calculate epsilon for each data point, and probability.
	for(z in 1:n) {
		eps[z] <- log(time[z]/tpred[z]) #epsilon: error for that point
		eps_back[z] <- log(time_back[z]/tpred[z]) #epsilon minus 1
		p[z] <- (1 - censored[z])*(phi((eps[z] + 0.5*sigma_e^2)/sigma_e) - phi((eps_back[z] + 0.5*sigma_e^2)/sigma_e)) + censored[z]*(1 - phi((eps[z] + 0.5*sigma_e^2)/sigma_e))
	}
	
	# Likleihood
	
	for(z in 1:n) {
		y[z] ~ dbern(max(10^-10, p[z]))
	}
	
} # end MTE model
# ###############################################################################
# # logLikelihood version for optim() in R
# ###############################################################################
logLikDevelop <- function(params, dat.list, mod){
	
	o <- names(params)
	
		
	tpred <- numeric(dat.list$n)
	
	#----------------------------------
	if(mod == "IT"){
		
		for(z in 1:dat.list$n) {
			tpred[z] <- params[which(o == paste("tau[", dat.list$temp[z], ",", dat.list$loc[z], "]", sep=""))]#*params[which(o== paste("v[", dat.list$incubator[z], "]", sep=""))] 
		}
	
	#----------------------------------
	} else if(mod == "SSL"){
		
		tau <- matrix(NA, nrow = dat.list$nT, ncol = dat.list$nL)
		for(i in 1:dat.list$nT){
			for(l in 1:dat.list$nL){
				tau[i,l] <- params[paste("tau0[", l, "]", sep="")] * exp(params[paste("E[", l, "]", sep="")]/(8.62*10^-5) * (1/(dat.list$T.obs[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(params[paste("El[", l, "]", sep="")]/(8.62*10^-5) * (1/(dat.list$T.obs[i] + 273.15) - 1/(params[paste("Tl[", l, "]", sep="")] + 273.15)))) 
			}}
		
		tpred <- numeric(dat.list$n)
		for(z in 1:dat.list$n) {
			tpred[z] <- tau[dat.list$temp[z], dat.list$loc[z]]#*params[which(o== paste("v[", dat.list$incubator[z], "]", sep=""))]
			}
	
	#----------------------------------
	} else if(mod == "SSUL"){
		
		tau <- matrix(NA, nrow = dat.list$nT, ncol = dat.list$nL)
		for(i in 1:dat.list$nT){
			for(l in 1:dat.list$nL){
				tau[i,l] <- params[paste("tau0[", l, "]", sep="")] * exp(params[paste("E[", l, "]", sep="")]/(8.62*10^-5) * (1/(dat.list$T.obs[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(params[paste("El[", l, "]", sep="")]/(8.62*10^-5) * (1/(dat.list$T.obs[i] + 273.15) - 1/(params[paste("Tl[", l, "]", sep="")] + 273.15))) + exp(params[paste("Eh[", l, "]", sep="")]/(8.62*10^-5) * (-1/(dat.list$T.obs[i] + 273.15) + 1/(params[paste("Th[", l, "]", sep="")] + 273.15)))) 
			}}
		
		tpred <- numeric(dat.list$n)
		for(z in 1:dat.list$n) {
			tpred[z] <- tau[dat.list$temp[z], dat.list$loc[z]]#*params[which(o== paste("v[", dat.list$incubator[z], "]", sep=""))]
		}
	
		#----------------------------------
	} else if(mod == "SSUL_1L"){
		
		tau <- matrix(NA, nrow = dat.list$nT, ncol = dat.list$nL)
		for(i in 1:dat.list$nT){
			tau[i,1] <- params["tau0"] * exp(params["E"]/(8.62*10^-5) * (1/(dat.list$T.obs[i] + 273.15) - 1/(15 + 273.15))) * (1 + exp(params["El"]/(8.62*10^-5) * (1/(dat.list$T.obs[i] + 273.15) - 1/(params["Tl"] + 273.15))) + exp(params["Eh"]/(8.62*10^-5) * (-1/(dat.list$T.obs[i] + 273.15) + 1/(params["Th"] + 273.15)))) 
			}
		
		tpred <- numeric(dat.list$n)
		for(z in 1:dat.list$n) {
			tpred[z] <- tau[dat.list$temp[z], dat.list$loc[z]]#*params[which(o== paste("v[", dat.list$incubator[z], "]", sep=""))]
		}
		
		#----------------------------------
	} else if(mod == "A"){
		
		tau <- matrix(NA, nrow = dat.list$nT, ncol = dat.list$nL)
		for(i in 1:dat.list$nT){
			for(l in 1:dat.list$nL){
				tau[i,l] <- params[paste("tau0[", l, "]", sep="")] * exp(params[paste("E[", l, "]", sep="")]/(8.62*10^-5) * (1/(dat.list$T.obs[i] + 273.15) - 1/(15 + 273.15)))
			}}
		
		tpred <- numeric(dat.list$n)
		for(z in 1:dat.list$n) {
			tpred[z] <- tau[dat.list$temp[z], dat.list$loc[z]]#*params[which(o== paste("v[", dat.list$incubator[z], "]", sep=""))]
		}
		
		#----------------------------------
	}
	
	# Calculate epsilon for each data point, and probability.
	eps <- numeric(dat.list$n)
	eps_back <- numeric(dat.list$n)
	p <- numeric(dat.list$n)
	for(z in 1:dat.list$n) {
		eps[z] <- log(dat.list$time[z]/tpred[z]) #epsilon: error for that point
		eps_back[z] <- log(dat.list$time_back[z]/tpred[z]) #epsilon minus 1
		
		p[z] <- (1 - dat.list$censored[z]) * (pnorm((eps[z] + 0.5*params['sigma_e']^2)/params['sigma_e']) - pnorm((eps_back[z] + 0.5*params['sigma_e']^2)/params['sigma_e'])) 
			+ dat.list$censored[z]*(1 - pnorm((eps[z] + 0.5*params['sigma_e']^2)/params['sigma_e']))
	}
	
	# Likleihood
	LLv <- dnorm(params[which(o=="v[1]"):which(o=="v[56]")], mean = 1, sd = params['sigma_v'], log=TRUE)
	LLd <- numeric(dat.list$n)
	for(z in 1:dat.list$n) {
		LLd[z] <- dbinom(1, size = 1, prob = max(p[z], 10^-10), log=TRUE)
	}

	return(list(NLL = -sum(LLd, LLv), -LLd, -LLv))
		}