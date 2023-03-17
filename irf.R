options(rgl.useNULL=TRUE)
library(iRF)
library(AUC)

source('./fileio.R')
source('./id_mapping.R')

get_important_features.irf <- function (data_file, response, ref_level, mask = NULL) {
  ds <- load.dataset(
    meta.file = './data/sample_sheet.csv', meta.sep = ',',
    data.file = data_file, data.sep = ','
  )
  
  if (!is.null(mask)) {
    ds$data.matrix <- ds$data.matrix[, colnames(ds$data.matrix) %in% mask]
  }
  
  X <- as.data.frame(ds$data.matrix)
  Y <- ds$meta.data[[response]]
  
  # prepare inputs
  X <- X[Y != 'Control',]
  Y <- as.factor(Y[Y != 'Control'])
  Y <- as.factor(as.numeric(relevel(Y, ref_level)) - 1)
  
  
  n <- dim(X)[1]
  p <- dim(X)[2]
  
  # split training-testing set
  train.id <- sample(1:n, size = 8 * round(n / 10))
  test.id <- setdiff(1:n, train.id)
  
  # fit iRF without iterations
  # sel.prob <- rep(1/p, p)
  # 
  # fit <- iRF(x = X[train.id,], 
  #            y = Y[train.id], 
  #            xtest = X[test.id,], 
  #            ytest = Y[test.id],
  #            n.iter = 5, 
  #            iter.return = 1:5,
  #            n.core = 4
  # )
  
  # plot ROC/AUC
  # rf <- fit$rf.list
  
  # plot(0:1, 0:1, type = 'l', lty = 2, xlab = 'FPR', ylab = 'TPR', main = 'ROC Curve')
  # for (iter in 1:5){
  #   cat(paste('iter = ', iter, ':: '))
  #   roc.info <- roc(rf[[iter]]$test$votes[,2], Y[test.id])
  #   lines(roc.info$fpr, roc.info$tpr, type = 'l', col = iter, lwd = 2)
  #   cat(paste('AUROC: ', round(100*auc(roc.info), 2), '%\n', sep = ''))
  # }
  # legend('bottomright', legend=paste('iter:', 1:iter), col=1:iter, lwd=2, bty='n')
  
  
  # plot feature weights
  # par(mfrow = c(1,5))
  # for (iter in 1:5){
  #   varImpPlot(rf[[iter]], n.var = 10, main = paste('Variable Importance (iter:', iter, ')')) 
  # }
  
  
  # fit iRF with iterations
  fit <- iRF(
    x = X[train.id,],
    y = Y[train.id],
    xtest = X[test.id,], 
    ytest = Y[test.id],
    n.iter = 10,
    iter.return = 1:10,
    n.bootstrap = 50,
    select.iter = T,
    n.core = 6
  )
  
  features_importance <- as.data.frame(fit$rf.list$importance)
  significant_features <- features_importance %>%
    arrange(desc(MeanDecreaseGini)) %>%
    filter(MeanDecreaseGini != 0)
  
  return(list(
    significant_features = significant_features,
    irf_fit = fit
  ))
}

plot_stability_scores <- function (irf_fit, filename) {
  toplot <- irf_fit$interaction$stability
  names(toplot) <- irf_fit$interaction$int
  toplot <- sort(toplot, decreasing = T)
  
  png(filename = filename)
  dotchart(rev(toplot[1:min(20, length(toplot))]),
           xlab = 'Stability Score',
           xlim = c(0, 1)
  )
  dev.off()
}


