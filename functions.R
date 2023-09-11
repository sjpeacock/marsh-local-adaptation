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

#------------------------------------------------------------------------------
# Function to predict parameter value for different MTE functions
# @param temp vector of temperatures for which parameter value will be returned
# @param params named vector of hypoerparameters associated with the MTE model
# @param fn the MTE function to be applied (A, SSU, SSUL, or SSL)
# @returns vector of parameter value (tau or mu) for each temperatured in temp
#------------------------------------------------------------------------------

predictMTE <- function(temp, params, fn="A"){
	if(length(which(names(params) == 'mu0')) == 0){
		a <- params['tau0']
		Efront <- 1
	} else {
		a <- params['mu0']
		Efront <- -1
	}
	if(fn=="A"){
		parEst <- a * exp(Efront * params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15)))
		
	} else if (fn=="SSU"){
		parEst <- a * exp(Efront * params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['Eh']/(8.62*10^-5) * (-1/(temp + 273.15) + 1/(params['Th'] + 273.15))))
		
	} else if (fn=="SSUL"){
		parEst <- a * exp(Efront * params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['El']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(params['Tl'] + 273.15))) + exp(params['Eh']/(8.62*10^-5) * (-1/(temp + 273.15) + 1/(params['Th'] + 273.15))))
		
	} else if (fn=="SSL"){
		parEst <- a * exp(Efront * params['E']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(15 + 273.15))) * (1 + exp(params['El']/(8.62*10^-5) * (1/(temp + 273.15) - 1/(params['Tl'] + 273.15))))
		
	} else {
		stop("Function not defined.")
	}
	return(parEst)
}

#------------------------------------------------------------------------------
# Function to predict R0 over temperature assuming SSUL for tau and SSU for mu
# @param temp vector 
# @param params named list (mu, tau) with each element containing a named list of the 
#        hypoerparameters associated with the MTE model for that parameter
# @returns vector of parameter value (tau or mu) for each temperatured in temp
#------------------------------------------------------------------------------

R0calc <- function(temp, params, DL.out = FALSE, rhoH = NULL){
	mu <- predictMTE(temp = temp, params = params$mu, fn = "SSU")
	tau <- predictMTE(temp = temp, params = params$tau, fn = "SSUL")
	DL <- exp(- mu * tau)
	
	if(DL.out == TRUE){
		return(DL)
	} else {
		R0 <- DL * rhoH / (mu + rhoH)
		return(R0)
	}
}

#------------------------------------------------------------------------------
# Function to calculate DL, the proportion of eggs surviving to infectivity,
# over temperature assuming SSUL for tau and SSU for mu
# @param temp vector 
# @param params named list (mu, tau) with each element containing a named list of the 
#        hypoerparameters associated with the MTE model for that parameter
# @returns vector of parameter value (tau or mu) for each temperatured in temp
#------------------------------------------------------------------------------

DLcalc <- function(temp, params){
	mu <- predictMTE(temp = temp, params = params$mu, fn = "SSU")
	tau <- predictMTE(temp = temp, params = params$tau, fn = "SSUL")
	DL <- exp(- mu * tau)
	return(DL)
}



