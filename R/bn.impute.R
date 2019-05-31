#' @import bnlearn

bn.impute <- function(exp.dat, models, genes=NULL){
  if (is.null(genes))
    genes <- names(models)
  
  predicted.values <- c()
  for (g in seq_len(length(genes))){
    sub.dat <- as.data.frame(exp.dat[,colnames(exp.dat) %in% nodes])
    predicted.values <- c(predicted.values, predict(training[g], node=genes[g], data=sub.dat, method='bayes-lw'))
  }
  return(res)
}

