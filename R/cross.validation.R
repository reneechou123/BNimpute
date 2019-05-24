#' @import parallel
#' @import bnlearn
#' @import Rgraphviz
#' @import WGCNA
#' @import dplyr
#' @import rlist

cross.validation <- function(ref.exp, modules, eig.exp, power, gene.exp.min=0.5, gene.exp.max=10, num.genes=3, seed=123, ..., replace.unidentifiable=TRUE){

  # validation function
  NRMSE <- function(true.data, imputed.data){
    if (sum(true.data != 0) == 0 & sum(imputed.data != 0) == 0)
      return(0)
    return(sqrt(sum((true.data - imputed.data) ^ 2) / sum((true.data - mean(true.data)) ^ 2)))
  }

  # preprocessing
  low.genes <- colnames(ref.exp[,apply(ref.exp, 2, function(x) sum(x > gene.exp.max) == 0 & sum(x < gene.exp.min) == 0)])
  eig.exp <- eig.exp[,colnames(eig.exp)!='MEgrey']
  modules[,1] <- paste0('ME', modules[,1])
  num.modules <- length(unique(modules[,1])) - 1 # minus grey module

  # select highly connected genes as model variables
  res <- lapply(seq_len(num.genes), function(x) chooseOneHubInEachModule(ref.exp, colorh=modules[,1], numGenes=100,
                                                                  omitColors='MEgrey', power=power))
  res <- as.data.frame(matrix(unlist(res), ncol=num.modules, byrow=TRUE, dimnames=list(NULL, names(res[[1]]))),
                       stringAsFactors=FALSE)

  set.seed(seed)
  error.values <- c()
  predicted.matrix <- matrix(, nrow=length(low.genes), ncol=0)
  rownames(predicted.matrix) <- low.genes
  for (n in seq_len(nrow(ref.exp))){
    
    # build models for low expression genes
    print(paste0('Cross validation iteration ', n, '...'))
    training.data <-ref.exp[-n,,drop=F]
    test.data <- ref.exp[n,,drop=F]
    predicted.values <- c()
    predicted.genes <- c()

    for (g in seq_len(length(low.genes))){
      gene <- low.genes[g]
      bn.dat <- cbind(ref.exp[,colnames(ref.exp)==gene,drop=FALSE], eig.exp[rownames(eig.exp) %in% rownames(ref.exp),])
      bn.dat[[gene]] <- scale(bn.dat[[gene]])
      if (sum(is.na(bn.dat)) > 0){
        next
      }

      # construct network
      structure.1 <- hc(bn.dat)
      if (!gene %in% structure.1$arcs){
        print(paste0('Gene ', gene, ' emitted during hill-climbing structure learning process.'))
        next
      }

      from <- c()
      to <- c()
      sub.res <- res
      sub.res[sub.res == gene] <- NA
      for (i in seq_len(nrow(structure.1$arcs))){
        f <- as.character(sub.res[[structure.1$arcs[i, 1]]]) # from
        c <- as.character(sub.res[[structure.1$arcs[i, 2]]]) # to
        if (length(f)==0){
          from <- c(from, rep(structure.1$arcs[i, 1], length(c)))
          to <- c(to, c)
        }
        else if(length(c)==0){
          from <- c(from, f)
          to <- c(to, rep(structure.1$arcs[i, 2], length(f)))
        }
        else{
          from <- c(from, f)
          to <- c(to, c)
        }
      }
      arc.set <- unique(matrix(c(from, to), nrow=length(from)))
      arc.set <- arc.set[arc.set[,1] != arc.set[,2],]
      arc.set <- arc.set[complete.cases(arc.set),]
      nodes <- as.character(unique(c(from, to)))
      nodes <- nodes[!is.na(nodes)]
      structure.2 <- empty.graph(nodes)
      arcs(structure.2) <- arc.set
      
      # train the model
      sub.exp <- as.data.frame(ref.exp[,colnames(ref.exp) %in% nodes])
      training <- bn.fit(structure.2, sub.exp, ..., replace.unidentifiable=replace.unidentifiable)
      
      # estimate missing vlaues
      sub.test <- as.data.frame(test.data[,colnames(test.data) %in% nodes])
      predicted.values <- c(predicted.values, predict(training, node=gene, data=sub.test, method='bayes-lw'))
      predicted.genes <- c(predicted.genes, gene)
    }

    true.data <- t(as.matrix(test.data[,colnames(test.data) %in% predicted.genes]))
    predicted.values <- as.matrix(predicted.values)
    rownames(predicted.values) <- predicted.genes
    nrmse <- NRMSE(true.data, imputed.data=predicted.values)
    error.values <- c(error.values, nrmse)
    predicted.matrix <- cbind(predicted.matrix, predicted.values[match(rownames(predicted.matrix), rownames(predicted.values))])
  }  

  colnames(predicted.matrix) <- row.names(ref.exp)
  res <- list(error.values, predicted.genes, predicted.matrix)
  names(res) <- c('error.values', 'predicted.genes', 'predicted.matrix')
  return(res)
}

