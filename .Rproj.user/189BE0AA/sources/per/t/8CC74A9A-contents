

rm(list=ls())


# Load the required R-libraries
# Standard filter package
library(mFilter)
# Multivariate time series: VARMA model for macro indicators: used here for simulation purposes only
library(MTS)
# HAC estimate of standard deviations in the presence of autocorrelation and heteroscedasticity
library(sandwich)
# Extended time series
library(xts)
# Library for Diebold-Mariano test of equal forecast performance
library(multDM)
# GARCH model: for improving regression estimates
library(fGarch)



# Load the relevant M-SSA functionalities
# M-SSA functions
source(paste(getwd(),"/R/functions_MSSA.r",sep=""))
# Load signal extraction functions used for JBCY paper (relies on mFilter)
source(paste(getwd(),"/R/HP_JBCY_functions.r",sep=""))
# Utility functions for M-SSA, see tutorial 
source(paste(getwd(),"/R/M_SSA_utility_functions.r",sep=""))
# Set of performance metrics and tests of unequal predictability
#source(paste(getwd(),"/R/performance_statistics_functions.r",sep=""))





 
path.main<-getwd()

path.pgm<-paste(path.main,"/R/",sep="")
path.out<-paste(path.main,"/Latex/",sep="")
path.sweave<-paste(path.main,"/Sweave/",sep="")
path.data<-paste(path.main,"/Data/",sep="")
# Savd results from an empirical analysis of S&P500
path.result<-paste(path.main,"/Results/",sep="")
fig_size<-4

#------------------------------

recompute_results<-T

paper<-"BIP_predictor"

script <- paste(path.sweave,paper,sep="")
## enforce par(ask=FALSE)
options(device.ask.default=FALSE)
## create a LaTeX file
Sweave(script,output=paste(path.out,paper,".tex",sep=""))

#Stangle(script)

