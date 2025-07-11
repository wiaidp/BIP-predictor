
# Created 25-05-2025 by Simon van Norden
# This code just loads the BIP data set.
# It is based on snippets from BIP_predictor.Rnw Section 1 written by Marc Wildi
#
GetBIPdata <- function(main_path = NULL, DataFrame = F) {
  # if main_path omitted, it defaults to ~/GitHub/BIP-predictor
  path.data <- ifelse( is.null(main_path),
                       "C:/Users/svn/Documents/GitHub/BIP-predictor/Data/",
                       main_path)
  
  #------------------------------
  
  # Load data
  # Transformed indicators: differences, trimmed, standardized
  load(file=paste0(path.data,"macro"))
  
  select_vec_multi<-c("BIP","ip","ifo_c","ESI","spr_10y_3m")
  x_mat<-data[,select_vec_multi] 
  rownames(x_mat)<-rownames(data)
  n<-dim(x_mat)[2]
  # Number of observations
  len<-dim(x_mat)[1]
  
  return(x_mat)
}


mat_lag <- function(X, p) {
  n <- nrow(X)
  k <- ncol(X)
  vnames <- colnames(X, do.NULL = F)
  # Initialize lagged matrix
  lagged_matrix <- matrix(NA, nrow = n - p, ncol = k * p)
  colnames(lagged_matrix) <- paste(rep(vnames,p), rep(1:p, each = k), sep = 'L')
  rownames(lagged_matrix) <- rownames(X[-(1:p), ])
  
  for (j in 1:p) {
    lagged_matrix[, ((j-1)*k + 1):(j*k)] <- X[(p-j+1):(n - j), ]
  }
  
  return(lagged_matrix)
  
}


enetVAR <- function(X, p = 1, nIRF = 20, nsteps = NULL, mode = 'step', 
                    quiet = F, cpalette = "ggplot2", lwd = 3, cex = 1){
  
  nvar <- ncol(X) 
  Xnames <- colnames(X)
  MyX <- mat_lag(X, p)  # Include 1:p lags of all variables 
  MyY <- X[-(1:p),]     # Drop the first p lags of all variables
  
  # Clear space to hold the enet objects
  MyENets <- vector("list", nvar) # list of enet() results
  names(MyENets) <- Xnames
  
  # Create list to hold returned objects
  out_list <- vector("list", 4)   # list of things to return
  names(out_list) <- c("Phi", "IRF", "Resids", "Sigma")
  out_list$Phi <- matrix(NA, nrow = nvar, ncol = p*nvar)
  rownames(out_list$Phi) <- Xnames
  colnames(out_list$Phi) <- paste(rep(Xnames,p),
                                  rep(1:p, each = nvar),
                                  sep = '_L')
  
  # Clear space to hold the fitted values and residuals
  ENetFits <- matrix(NA, nrow = nrow(X)-p, ncol = nvar)
  colnames(ENetFits) <- Xnames
  
  # default value for degree of regularization should be LS
  if (is.null(nsteps)) nsteps = 1 + ncol(MyX)
  # Estimate the ElasticNets for each equation of the VAR
  for (j in 1:nvar) {
    # Estimate the enet 
    MyENets[[j]] <- enet(MyX, MyY[,j])
    # Store the regularized coefficients 
    out_list$Phi[j,] <- MyENets[[j]]$beta.pure[nsteps,]
  }
  if (!quiet) {
    print("Estimated VAR Coefficients:")
    print(out_list$Phi)
  }
  
  # Calculate the IRFs for the *first* series only! 
  # To change this, change  $psi[1,]   to  $psi  on the next line 
  #   and redimension the results
  out_list$IRF <- matrix(VARpsi(out_list$Phi, lag = nIRF)$psi[1,], ncol = nvar, byrow = T)
  colnames(out_list$IRF) <- Xnames   # labels source of shocks
  rownames(out_list$IRF) <- seq.int(from = 0, to = nIRF) # labels lags     
  if (!quiet) {
    print("Estimated IRF for 1st Variable:")
    print(t(out_list$IRF))
    
    # save the old colour palette, just in case
    old_colours <- palette()
    if (!identical(old_colours, cpalette)) palette(cpalette)
    plot(0:nIRF, out_list$IRF[,1], type = "l",
         main = paste0("Regularized VAR(", p, "): ElasticNet"),
         ylab = paste("IRF for", Xnames[1]),
         xlab = "Lag",
         lwd = 3, col = 1)
    for (j in 2:nvar)
      lines(0:nIRF, out_list$IRF[,j], lwd = 3, col = j)
    legend("topright", 
           legend = Xnames, 
           col = 1:nvar, 
           lwd = 3,
           cex = cex,
           title = "Shocks to:")
    # restore the old colour palette
    palette(old_colours)
  }
  
  # Calculate the residuals and estimate Sigma
  for (j in 1:nvar) {
    # Retrive the fitted values for each equation
    ENetFits[,j] <- predict.enet(MyENets[[j]], 
                                 newx = MyX, 
                                 s = nsteps,
                                 type = "fit", 
                                 mode = mode)$fit
  }
  # Calculate the residuals
  out_list$Resids <- MyY - ENetFits
  # Estimate Sigma as the sample VCV of the residuals
  out_list$Sigma <- cov(out_list$Resids)
  
  if (!quiet) {
    print(paste('The VCV matrix of the estimated VAR residuals is:\n'))
    print( out_list$Sigma)
  }
  
  return(out_list)
  
}


library(elasticnet)


main_path<-paste(getwd(),"/Data/",sep="")
MyData <-  GetBIPdata(main_path)
tail(MyData)


p <- 4  # number of lags to include
MyX <- mat_lag(MyData[,-1], p)  # Include 1:p lags of all variables but the 1st
MyY <- MyData[-(1:p),1]         # Drop the first p lags of the 1st variable

# Estimate the enet for ALL POSSIBLE values
# Includes intercept, all variables normalized before estimate
enet_BIP <- enet(MyX, MyY)
# plot the results
plot(enet_BIP, use.color = T, lwd=3)

print(enet_BIP)

# show the estimated coefficients at each step in the path.
K <- ncol(MyX)
predict(enet_BIP, type = "coefficients", s = 2:(K+1))


# What's the order of the VAR we're modeling?
p<-3
p <- 2
# What's the depth of the ENet we want to use?
nsteps <- 6
nsteps <- 5
# How many variables in our VAR?
nvar <- ncol(MyData) 

# Clear space to hold the enet results
MyENets <- vector("list", nvar)
names(MyENets) <- colnames(MyData)

# Clear space for the Phi matrix
VARcoefs <- matrix(NA, nrow = nvar, ncol = p*nvar)

MyX <- mat_lag(MyData, p)  # Include 1:p lags of all variables 
MyY <- MyData[-(1:p),]         # Drop the first p lags of all variables
for (j in 1:nvar) {
  
  # Estimate the enet 
  MyENets[[j]] <- enet(MyX, MyY[,j])
  # Store the regularized coefficients 
  VARcoefs[j,] <- predict.enet(MyENets[[j]], type = "coefficients", s = nsteps,
                               mode = "step")$coefficients
}

print(t(VARcoefs))



# What do the IRFs look like for BIP?
IRFenet <- matrix(VARpsi(VARcoefs, lag = 20)$psi[1,], ncol = 5, byrow = T)
colnames(IRFenet) <- colnames(MyData)   # labels source of shocks
rownames(IRFenet) <- seq.int(from = 0, to = 20) # labels lags     

MyColours <- c("blue", "darkgreen", "red3", "gray80", "lightgreen", "orange")
plot(0:20, IRFenet[,1], type = "l",
     main = paste0("Regularized VAR(", p, "): ElasticNet"),
     ylab = "IRF for GDP",
     xlab = "Lag",
     lwd = 3, col = MyColours[1])
lines(0:20, IRFenet[,2], lwd = 3, col = MyColours[2])
lines(0:20, IRFenet[,3], lwd = 3, col = MyColours[3])
lines(0:20, IRFenet[,4], lwd = 3, col = MyColours[4])
lines(0:20, IRFenet[,5], lwd = 3, col = MyColours[5])
legend("topright", 
       legend = colnames(MyData), 
       col = MyColours[1:5], 
       lwd = 3,
       title = "Shocks to:")


nIRF = 20
quiet<-F

list_obj<-enetVAR(MyData, p, nIRF, nsteps , mode = 'step',quiet )

list_obj$Phi
list_obj$IRF
list_obj$Resids
list_obj$Sigma

  