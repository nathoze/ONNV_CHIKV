
library(questionr)
library(rstan)
library(rmutil)
library(MASS)
library(reshape2)
library(binom)
library(dplyr)

# Load data and functions ----
source('loadAll.R')

## MCMC  ----

source('launch_MCMC.R')

## extract the chains 
Chains= rstan::extract(fit)


# Estimate the DIC of the model
source('DIC_stan_model.R')


## posterior distributions of parameters 
source('parameter_estimates.R')
 
# seroprevalence using the classical method
source('ClassicalSeroprevalence.R')


## Figure 1B
source('Figure1.R')


# # Table 4 
source('Compare_simulations_observations.R')
C =  ComparePosteriorDistObs(data_csv, Chains, regions=1)

V=array(data=0,dim = c(5,5))

for( i in seq(0,4)){
  for(j in 1:5){
    P=paste(C[2*i+1,j],C[2*i+2,j])
    V[i+1,j] =P 
  }
}

