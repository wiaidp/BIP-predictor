# ROCplots : a function for comparing AUCs and ROC plots
#   by Simon van Norden, HEC Montreal
#   19-05-2025
#
# SYNTAX
#   
#   auc <- ROCplots(df, showROC = TRUE, smoothROC = FALSE, showLegend = TRUE)
#
# where
#   df  is a data.frame of series to analyse. See below for format.
#   showROC  controls creation of the plot comparing ROC curves
#   smooth   controls whether ROC curves should be smoothed
#   showLegend controls where a legend is added to the plot of ROC curves
#
# VALUE
#
#   auc   a 1-column data.frame containing the AUC for each predictor in df.
#
# NOTES
#
# The data.frame df should contain
# - the target series in column 1, consisting only of 0s and 1s
# - an arbitrary number of additional numeric series in the remaining columns
# All series should be the same length and contain no NAs

inst_pack<-rownames(installed.packages())
if (!"pROC"%in%inst_pack)
install.packages("pROC")

library(pROC)

ROCplots <- function(df, showROC = T, smoothROC = FALSE, showLegend = TRUE) {
  df_names <- names(df) # series names
  k <- length(df_names) # 1 + number of predictors to compare
  print(paste("Target variable is", df_names[1]))
  # Create table to hold AUCs
  aurocs <- as.data.frame(matrix(NA, nrow = k-1))
  row.names(aurocs) <- df_names[-1]
  names(aurocs) <- "AUC"
  
  # make the first plot and store the first AUC
  ROC_OBJ <- roc_(df, response = df_names[1], 
                  predictor = df_names[2],
                  quiet = T, plot = showROC,
                  smooth = smoothROC)
  aurocs[1,1] <- ROC_OBJ$auc
  
  if (k>2) {
  # Calculate and plot for the rest of the series
    for (j in 3:k) {
      ROC_OBJ <- roc_(df, response = df_names[1], 
                      predictor = df_names[j], 
                      quiet = T, plot = F,
                      smooth = smoothROC)
      aurocs[j-1,1] <- ROC_OBJ$auc
      # Plot the ROC?
      if (showROC) lines(ROC_OBJ$specificities, 
                         ROC_OBJ$sensitivities, 
                         col = j-2, lwd = 2)
    }
  }
  
  if (showROC){
    if (k > 5) {
      roc_cex = 0.75
      roc_ncol = 2
    } else {
      roc_cex = 1
      roc_ncol = 1
    }
    if (showLegend) legend("bottomright", 
                           legend = df_names[-1], 
                           col = 1:(k-1), 
                           lwd = 3,
                           cex = roc_cex,
                           ncol = roc_ncol)
  }
  return(aurocs)
}
