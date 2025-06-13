


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


