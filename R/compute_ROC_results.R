# ROC curve sample code with artificial data

sd <- seq(from = 1, to = 0.1, by = -0.1) # Std Dev's of Errors
k <- length(sd) # No. of Series
n <- 100 # No. of Observations

# Generate target and k series of length n
set.seed(20250519)
mytarget <- matrix(rnorm(n), nrow = n)
mypreds  <- matrix(rnorm(n*k, mean = 0, sd = sd), 
                   nrow = n, byrow = T)
# Predictors = Target + noise
mypreds <- matrix(mytarget, ncol = k, nrow = n) + mypreds

plot(mytarget,mypreds[,1], main = "Demo Data",
     xlab = "Target", ylab = "Predictors") 
for (j in 2:10) points(mytarget, mypreds[,j], col = j)

Mydata <- data.frame(cbind(mytarget,mypreds))
names(Mydata) <- c("Target", paste0("P_",1:10))
# convert Target >0 to 1, rest is 0
Mydata$Target <- ifelse(Mydata$Target > 0, 1, 0 )



ROCplots(Mydata, smoothROC = T)

########################################################################################################################
# Alternative BVAR(4) model specification

# Try a BVAR(4)
V0_BVAR <- diag(rep(1,n))  # prior for Sigma_a matrix - K x K

nlags <- 24
lambda_BVAR <- 10  # scaling factor for precision matrix
p_BVAR <- 4
C_BVAR <- lambda_BVAR*diag(rep(1,k*p_BVAR+1))
BVAR4_obj <- BVAR(x_mat_wc,
                  p = p_BVAR, 
                  C = C_BVAR,
                  V0 = V0_BVAR)

BVAR4_MA <- VARpsi(Phi = BVAR4_obj$Phi, lag = nlags)
IRF_BIP_BVAR4 <- as.data.frame(t(matrix(BVAR4_MA$psi[1,], ncol = 1+nlags)))
names(IRF_BIP_BVAR4) <- labels(x_mat)[[2]]   # labels source of shocks
rownames(IRF_BIP_BVAR4) <- seq.int(from = 0, to = nlags) # labels lags

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


#########################################################################################################
# Apply ROC to direct forecast, direct HP forecast, M-MSE and M-SSA
# Data: 
#   -In-sample vs. out-of-sample
#   -Same sample (account for different NAs of various designs)
# Target: 
#   -GDP vs. HP-GDP
#   -forward-shift against h

# 1. Full sample without pandemic; all explanatories

# 1.1. Compute direct forecasts for h in h_vec

select_direct_indicator<-select_vec_multi
select_direct_indicator<-c("ifo_c","ESI")
# Select data matrix
data_roc<-x_mat_wc
# Note: too complex designs (too many indicators) lead to overfitting and thus worse out-of-sample performances
# To illustrate the direct predictor consider the following example of a h-step ahead direct forecast:
direct_forecast_mat<-NULL
for (h in 0:max(h_vec))
{
  # Shift GDP forward by publication lag+forecast horizon
  forward_shifted_GDP<-c(data_roc[(1+lag_vec[1]+h):nrow(data_roc),"BIP"],rep(NA,h+lag_vec[1]))
  # Regress selected indicators on forward-shifted GDP: start at t=L (same sample as HP-C)
  lm_obj<-lm(forward_shifted_GDP[L:length(forward_shifted_GDP)]~data_roc[L:nrow(data_roc),select_direct_indicator])
  sum_obj<-sum_obj_unfiltered<-summary(lm_obj)
  # Compute the predictor: one can rely on the generic R-function predict or compute the predictor manually
  direct_forecast<-lm_obj$coef[1]+data_roc[,select_direct_indicator]%*%lm_obj$coef[2:(length(select_direct_indicator)+1)]
  # Note that this is a full-sample predictor (no out-of-sample span)
  
  # We can now plot target and direct forecast: for h>2 the predictor comes close to a flat line centered at zero
  direct_forecast_mat<-cbind(direct_forecast_mat,direct_forecast)
}
colnames(direct_forecast_mat)<-paste("h=",h_vec,sep="")

#-----------
# 1.2. Direct HP forecast
direct_hp_forecast_mat<-forward_shifted_GDP_mat<-NULL
for (h in 0:max(h_vec))
{
  # Shift GDP forward by publication lag+forecast horizon
  forward_shifted_GDP<-c(data_roc[(1+lag_vec[1]+h):nrow(data_roc),"BIP"],rep(NA,h+lag_vec[1]))
  forward_shifted_GDP_mat<-cbind(forward_shifted_GDP_mat,forward_shifted_GDP)
  hp_mat<-NULL
  for (i in 1:ncol(data_roc))
   hp_mat<-cbind(hp_mat,filter(data_roc[,i],hp_c,side=1))
  colnames(hp_mat)<-select_vec_multi
  lm_obj<-lm(forward_shifted_GDP~hp_mat[,select_direct_indicator])
  # Compute the predictor: one can rely on the generic R-function predict or compute the predictor manually
  direct_hp_forecast_mat<-cbind(direct_hp_forecast_mat,lm_obj$coef[1]+hp_mat[,select_direct_indicator]%*%lm_obj$coef[2:(length(select_direct_indicator)+1)])
}
colnames(direct_hp_forecast_mat)<-colnames(forward_shifted_GDP_mat)<-paste("h=",h_vec,sep="")

#--------------
# 1.3. M-SSA and M-MSE
final_mssa_array<-final_mssa_indicator_obj$mssa_array
mssa<-final_mssa_array["BIP",,]
final_mmse_array<-final_mssa_indicator_obj$mmse_array
mmse<-final_mmse_array["BIP",,]

mssa_mat<-mmse_mat<-NULL
for (h in 0:max(h_vec))
{
  # Shift GDP forward by publication lag+forecast horizon
  forward_shifted_GDP<-c(data_roc[(1+lag_vec[1]+h):nrow(data_roc),"BIP"],rep(NA,h+lag_vec[1]))
  forward_shifted_GDP_mat<-cbind(forward_shifted_GDP_mat,forward_shifted_GDP)
  lm_obj<-lm(forward_shifted_GDP~t(final_mssa_array[select_direct_indicator,,h+1]))
# Compute the predictor: one can rely on the generic R-function predict or compute the predictor manually
  mssa_mat<-cbind(mssa_mat,lm_obj$coef[1]+t(final_mssa_array[select_direct_indicator,,h+1])%*%lm_obj$coef[2:(length(select_direct_indicator)+1)])
  lm_obj<-lm(forward_shifted_GDP~t(final_mmse_array[select_direct_indicator,,h+1]))
  mmse_mat<-cbind(mmse_mat,lm_obj$coef[1]+t(final_mmse_array[select_direct_indicator,,h+1])%*%lm_obj$coef[2:(length(select_direct_indicator)+1)])
}
colnames(mssa_mat)<-colnames(mmse_mat)<-paste("h=",h_vec,sep="")


#-------------
# 1.4. Targets: forward-shifted GDP and HP-GDP
forward_shifted_GDP_mat<-forward_shifted_GDP_mat
forward_shifted_HP_GDP_mat<-NULL
for (h in 0:max(h_vec))
{
  # Shift GDP forward by publication lag+forecast horizon
  forward_shifted_GDP<-c(data_roc[(1+lag_vec[1]+h):nrow(data_roc),"BIP"],rep(NA,h+lag_vec[1]))
  forward_shifted_HP_GDP_mat<-cbind(forward_shifted_HP_GDP_mat,filter(forward_shifted_GDP,hp_two,side=2))
}
colnames(forward_shifted_HP_GDP_mat)<-paste("shift",h_vec,sep="")

#-------------------------
# Apply ROC

shift_vec<-0:5
hh_vec<-0:6
AUC_array<-array(dim=c(length(shift_vec),length(hh_vec),4))
dimnames(AUC_array)<-list(paste("shift=",shift_vec,sep=""),paste("h=",hh_vec,sep=""),c("Direct forecast","Direct HP forecast","M-MSE","M-SSA"))
par(mfrow=c(length(shift_vec),length(hh_vec)))
for (i in 1:length(shift_vec))
{ 
  for (j in 1:length(hh_vec))
  {
    shift<-shift_vec[i]
    h<-hh_vec[j]
    
# Select target
    target<-as.integer(forward_shifted_HP_GDP_mat[,shift+1]>0)
    target<-as.integer(forward_shifted_GDP_mat[,shift+1]>0)
    
    ROC_data<-cbind(target,direct_forecast_mat[,h+1],direct_hp_forecast_mat[,h+1],mmse_mat[,h+1],mssa_mat[,h+1])
    rownames(ROC_data)<-rownames(data_roc)
    colnames(ROC_data)<-c("Target","Direct forecast","Direct HP forecast","M-MSE","M-SSA")
    ROC_data<-as.data.frame(na.exclude(ROC_data))
    AUC<-ROCplots(ROC_data, smoothROC = T)
    AUC_array[i,j,]<-unlist(AUC)
  }
}

AUC_array["shift=0",,]
AUC_array["shift=1",,]
AUC_array["shift=2",,]
AUC_array["shift=3",,]
AUC_array["shift=4",,]
AUC_array["shift=5",,]

#ts.plot(scale(ROC_data),col=c("grey","black","red","green","blue"))



