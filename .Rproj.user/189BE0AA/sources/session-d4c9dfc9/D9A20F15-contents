# Created 16-05-2025 by Simon van Norden
# This code is based on snippets from BIP_predictor.Rnw Sections 1 & 2
#   written by Marc Wildi
#
# It explores robustness of the IRFs to alternative estimated models using MTS


rm(list=ls())

# inst_pack<-rownames(installed.packages())
# if (!"mFilter"%in%inst_pack)
#   install.packages("mFilter")
# if (!"xts"%in%inst_pack)
#   install.packages("xts")
# if (!"MTS"%in%inst_pack)
#   install.packages("MTS")
# if (!"sandwich"%in%inst_pack)
#   install.packages("sandwich")
# if (!"multDM"%in%inst_pack)
#   install.packages("multDM")
# if (!"fGarch"%in%inst_pack)
#   install.packages("fGarch")
# if (!"xtable"%in%inst_pack)
#   install.packages("xtable")


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
# Library for tables
library(xtable)



# # Load the relevant M-SSA functionalities
# # M-SSA functions
# source(paste(getwd(),"/R/functions_MSSA.r",sep=""))
# # Load signal extraction functions used for JBCY paper (relies on mFilter)
# source(paste(getwd(),"/R/HP_JBCY_functions.r",sep=""))
# # Utility functions for M-SSA, see tutorial 
# source(paste(getwd(),"/R/M_SSA_utility_functions.r",sep=""))
# # Set of performance metrics and tests of unequal predictability
# #source(paste(getwd(),"/R/performance_statistics_functions.r",sep=""))

# main working directory should be ~/GitHub/BIP-predictor
path.main<-getwd()

path.pgm<-paste(path.main,"/R/",sep="")
path.out<-paste(path.main,"/Latex/",sep="")
path.sweave<-paste(path.main,"/Sweave/",sep="")
path.data<-paste(path.main,"/Data/",sep="")
# Savd results from an empirical analysis of S&P500
path.result<-paste(path.main,"/Results/",sep="")
#fig_size<-4

#------------------------------

recompute_results<-F

# Load data
data_file_name<-c("Data_HWI_2025_02.csv","gdp_2025_02.csv")
# Original (un-transformed) indicators
data_quarterly<-read.csv(paste(getwd(),"/Data/",data_file_name[2],sep=""))
BIP_original<-data_quarterly[,"BIP"]
# Transformed indicators: differences, trimmed, standardized
load(file=paste(getwd(),"\\Data\\macro",sep=""))

select_vec_multi<-c("BIP","ip","ifo_c","ESI","spr_10y_3m")
x_mat<-data[,select_vec_multi] 
rownames(x_mat)<-rownames(data)
n<-dim(x_mat)[2]
# Number of observations
len<-dim(x_mat)[1]

# Settings for M-SSA
ht_mssa_vec<-c(6.380160,  6.738270,   7.232453,   7.225927,   7.033768)
f_excess<-c(5,rep(4,length(select_vec_multi)-1))
# In-sample span for VAR
date_to_fit<-"2008"
# Model orders
p<-1
q<-0
# Filter length
L<-31
# Publication lag: one quarter for BIP
lag_vec<-c(1,rep(0,ncol(x_mat)-1))

#------------------------------------------------------
# Section 1: data
# Plots section data

# Plot BIP and diff-log BIP
# file = "./Figures/data.pdf"
# pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)
par(mfrow=c(1,1))

# Plot the data
# The real-time BIP (red) is lagging the target (black) by lag_vec[1] quarters (publication lag)
mplot<-x_mat
colo<-c("black",rainbow(ncol(data)-1))
main_title<-"Data: standardized diff-log, Pandemic outliers trimmed"
plot(mplot[,1],main=main_title,axes=F,type="l",xlab="",ylab="",col=colo[1],lwd=1,ylim=c(min(na.exclude(mplot)),max(na.exclude(mplot))))
mtext(colnames(mplot)[1],col=colo[1],line=-1)
for (i in 1:ncol(mplot))
{
  lines(mplot[,i],col=colo[i],lwd=1,lty=1)
  mtext(colnames(mplot)[i],col=colo[i],line=-i)
}
abline(h=0)
axis(1,at=c(1,4*1:(nrow(mplot)/4)),labels=rownames(mplot)[c(1,4*1:(nrow(mplot)/4))])
axis(2)
box()
# dev.off()



# Plot BIP and diff-log BIP
# file = "./Figures/data_lags.pdf"
# pdf(file = paste(path.out,file,sep=""), paper = "special", width = 6, height = 6)

par(mfrow=c(1,2))
# Plot the data
# The real-time BIP (red) is lagging the target (black) by lag_vec[1] quarters (publication lag)
mplot<-x_mat[which(rownames(x_mat)>=2007&rownames(x_mat)<=2011),]
colo<-c("black",rainbow(ncol(data)-1))
main_title<-"Financial Crisis"
plot(mplot[,1],main=main_title,axes=F,type="l",xlab="",ylab="",col=colo[1],lwd=1,ylim=c(min(na.exclude(mplot)),max(na.exclude(mplot))))
mtext(colnames(mplot)[1],col=colo[1],line=-1)
for (i in 1:ncol(mplot))
{
  lines(mplot[,i],col=colo[i],lwd=1,lty=1)
  mtext(colnames(mplot)[i],col=colo[i],line=-i)
}
abline(h=0)
axis(1,at=c(1,4*1:(nrow(mplot)/4)),labels=rownames(mplot)[c(1,4*1:(nrow(mplot)/4))])
axis(2)
box()
mplot<-x_mat[which(rownames(x_mat)>=2019&rownames(x_mat)<=2022),]
colo<-c("black",rainbow(ncol(data)-1))
main_title<-"Pandemic"
plot(mplot[,1],main=main_title,axes=F,type="l",xlab="",ylab="",col=colo[1],lwd=1,ylim=c(min(na.exclude(mplot)),max(na.exclude(mplot))))
mtext(colnames(mplot)[1],col=colo[1],line=-1)
for (i in 1:ncol(mplot))
{
  lines(mplot[,i],col=colo[i],lwd=1,lty=1)
  mtext(colnames(mplot)[i],col=colo[i],line=-i)
}
abline(h=0)
axis(1,at=c(1,4*1:(nrow(mplot)/4)),labels=rownames(mplot)[c(1,4*1:(nrow(mplot)/4))])
axis(2)
box()

# dev.off()

# Data without pandemic
x_mat_wc<-x_mat[which(rownames(x_mat)<2020|rownames(x_mat)>2021),]


#--------------------------------------------------------
# Section 2: dependence
#--------------------------------------------------------

# Set the data set, estimation & refinement parameters
# data_fit <- x_mat     # Full Sample
data_fit <- x_mat_wc  # Pre-Pandemic

# Set threshold higher to select simpler models
threshold <- 2.5

# Calculate IRFs for how many lags?
nlags <- 24

#--------------------------------------------------------
# Section 2.1: Examine refined VARMAs
#--------------------------------------------------------

# Search over VARMA[p,q] and compare AIC, BIC
p.max <- 6
q.max <- 1
VARMA.results <- as.data.frame(matrix(NA, nrow = p.max*(q.max+1), ncol = 6))
colnames(VARMA.results) <- c('p','q','Thresh','Nparam','aic','bic')
for (j in 1:p.max) {
  for (k in 0:q.max) {
    VARMA_obj<-VARMA(data_fit,p=j,q=k, prelim = T, details = F, thres = threshold)
    VARMA.results$p[(j-1)*(q.max+1) + (k + 1)] <- VARMA_obj$ARorder
    VARMA.results$q[(j-1)*(q.max+1) + (k + 1)] <- VARMA_obj$MAorder
    VARMA.results$Thresh[(j-1)*(q.max+1) + (k + 1)] <- threshold
    VARMA.results$Nparam[(j-1)*(q.max+1) + (k + 1)] <- sum(VARMA_obj$coef != 0)
    VARMA.results$aic[(j-1)*(q.max+1) + (k + 1)] <- VARMA_obj$aic
    VARMA.results$bic[(j-1)*(q.max+1) + (k + 1)] <- VARMA_obj$bic
  }
}
VARMA.results
# The results seem to prefer refined VARs over VARMAs
# AIC likes a 13 parameter VAR(4), BIC prefers a 7 parameter VAR(3)

# Let's recover the IRFs
VARMA_30<-VARMA(data_fit,p=3,q=0, prelim = T, details = F, thres = threshold)
VARMAirf_3 <- VARMAirf(Phi = VARMA_30$Phi, Theta = VARMA_30$Theta,
                       Sigma = VARMA_30$Sigma, lag = nlags, orth = F)

VARMA_40<-VARMA(data_fit,p=4,q=0, prelim = T, details = F, thres = threshold)
VARMAirf_4 <- VARMAirf(Phi = VARMA_40$Phi, Theta = VARMA_40$Theta,
                       Sigma = VARMA_40$Sigma, lag = nlags, orth = F)

# The following recovers the IRF for BIP (because BIP is row 1 in $psi)
IRF_BIP_VARMA3 <- as.data.frame(t(matrix(VARMAirf_3$psi[1,], ncol = 1+nlags)))
names(IRF_BIP_VARMA3) <- labels(x_mat)[[2]]   # labels source of shocks
rownames(IRF_BIP_VARMA3) <- seq.int(from = 0, to = nlags) # labels lags

IRF_BIP_VARMA4 <- as.data.frame(t(matrix(VARMAirf_4$psi[1,], ncol = 1+nlags)))
names(IRF_BIP_VARMA4) <- labels(x_mat)[[2]]   # labels source of shocks
rownames(IRF_BIP_VARMA4) <- seq.int(from = 0, to = nlags) # labels lags

#--------------------------------------------------------
# Section 2.2: Examine unrestricted VARs
#--------------------------------------------------------

# Uses same set of observations to compare all models
VARorder(data_fit, maxp = 6)
VARorder(data_fit, maxp = 8)
VARorder(data_fit, maxp = 12)
# In all three cases, LRstat and AIC want 4 lags but BIC and HQ prefer 1.

# Estimate the VAR parameters without restrictions
VAR1_obj <- VAR(x = data_fit, p = 1)
VAR4_obj <- VAR(x = data_fit, p = 4)

# Let's check out the MA coefficients implied by the two sets of estimates
# The IRFs are stored in 5x5 blocks of coefficients for lags 0 through 20
nlags <- 24
VAR1_MA <- VARpsi(Phi = VAR1_obj$Phi, lag = nlags)
VAR4_MA <- VARpsi(Phi = VAR4_obj$Phi, lag = nlags)

# The following recovers the IRF for BIP (because BIP is row 1 in $psi)
IRF_BIP_VAR1 <- as.data.frame(t(matrix(VAR1_MA$psi[1,], ncol = 1+nlags)))
names(IRF_BIP_VAR1) <- labels(x_mat)[[2]]   # labels source of shocks
rownames(IRF_BIP_VAR1) <- seq.int(from = 0, to = nlags) # labels lags

IRF_BIP_VAR4 <- as.data.frame(t(matrix(VAR4_MA$psi[1,], ncol = 1+nlags)))
names(IRF_BIP_VAR4) <- labels(x_mat)[[2]]   # labels source of shocks
rownames(IRF_BIP_VAR4) <- seq.int(from = 0, to = nlags) # labels lags

#--------------------------------------------------------
# Section 2.3: Estimate BVARs
#--------------------------------------------------------

# parameters for BVAR
k <- 5   # dimension of system
V0_BVAR <- diag(rep(1,k))  # prior for Sigma_a matrix - K x K
lambda_BVAR <- 0.001  # scaling factor for precision matrix
# precision matrix of size (kp+1) x (kp+1)  if there's a constant
#   low values <=> low precision

# Try a BVAR(3)
p_BVAR <- 3
C_BVAR <- lambda_BVAR*diag(rep(1,k*p_BVAR+1))
BVAR3_obj <- BVAR(data_fit,
                  p = p_BVAR, 
                  C = C_BVAR,
                  V0 = V0_BVAR)

BVAR3_MA <- VARpsi(Phi = BVAR3_obj$Phi, lag = nlags)
IRF_BIP_BVAR3 <- as.data.frame(t(matrix(BVAR3_MA$psi[1,], ncol = 1+nlags)))
names(IRF_BIP_BVAR3) <- labels(x_mat)[[2]]   # labels source of shocks
rownames(IRF_BIP_BVAR3) <- seq.int(from = 0, to = nlags) # labels lags

# Try a BVAR(4)
lambda_BVAR <- 10  # scaling factor for precision matrix
p_BVAR <- 4
C_BVAR <- lambda_BVAR*diag(rep(1,k*p_BVAR+1))
BVAR4_obj <- BVAR(data_fit,
                  p = p_BVAR, 
                  C = C_BVAR,
                  V0 = V0_BVAR)

BVAR4_MA <- VARpsi(Phi = BVAR4_obj$Phi, lag = nlags)
IRF_BIP_BVAR4 <- as.data.frame(t(matrix(BVAR4_MA$psi[1,], ncol = 1+nlags)))
names(IRF_BIP_BVAR4) <- labels(x_mat)[[2]]   # labels source of shocks
rownames(IRF_BIP_BVAR4) <- seq.int(from = 0, to = nlags) # labels lags


#--------------------------------------------------------
# Section 2.4: Compare the IRFs
#--------------------------------------------------------

MyColours <- c("blue", "darkgreen", "red3", "gray80", "lightgreen", "orange")
# plot the results
for (j in names(IRF_BIP_VAR1)) {
  plot(IRF_BIP_VAR1[[j]], type = 'l', col = MyColours[1], lwd = 3,
       ylim = c(-0.25, 1), 
       ylab = 'IRF of GDP', xlab = "Lags", main = paste("Shocks to", j))
  lines(IRF_BIP_VAR4[[j]], col = MyColours[2], lwd = 3)
  lines(IRF_BIP_VARMA3[[j]], col = MyColours[3], lwd = 3)
  lines(IRF_BIP_VARMA4[[j]], col = MyColours[4], lwd = 2)
  lines(IRF_BIP_BVAR3[[j]], col = MyColours[5], lwd = 3)
  lines(IRF_BIP_BVAR4[[j]], col = MyColours[6], lwd = 3)
  
  abline(h = 0, lty = 3)
  legend("topright", 
         legend = c("VAR(1)", "VAR(4)", "Reg. VAR(3)", "Reg. VAR(4)", 
                    "BVAR(3)", "BVAR(4)"), 
         col = MyColours, 
         lwd = 3)
}


