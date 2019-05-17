#' @import parallel
#' @import bnlearn
#' @import Rgraphviz
#' @import WGCNA
#' @import dplyr
#' @import rlist

build.model <- function(ref.exp, modules, eig.exp, power, gene.exp.min=0.5, gene.exp.max=10, num.genes=3, seed=123, ...){
  
  # preprocessing
  low.genes <- colnames(ref.exp[,apply(ref.exp, 2, function(x) sum(x > gene.exp.max) == 0 & sum(x < gene.exp.min) == 0)])
  eig.exp <- eig.exp[,colnames(eig.exp)!='MEgrey']
  modules[,1] <- paste0('ME', modules[,1])
  num.modules <- length(unique(modules[,1])) - 1 # minus grey module

  # select highly connected genes as model variables
  res <- lapply(1:num.genes, function(x) chooseOneHubInEachModule(ref.exp, colorh=modules[,1], numGenes=100, 
								  omitColors='MEgrey', power=power))
  res <- as.data.frame(matrix(unlist(res), ncol=num.modules, byrow=TRUE, dimnames=list(NULL, names(res[[1]]))), 
		       stringAsFactors=FALSE)

  # build models for low expression genes
  set.seed(123)
  models <- list()
  genes <- c()
  for (g in 1:length(low.genes)){
    gene <- low.genes[g]
    bn.dat <- cbind(ref.exp[,colnames(ref.exp)==gene,drop=FALSE], eig.exp[rownames(eig.exp) %in% rownames(ref.exp),])
    bn.dat[[gene]] <- scale(bn.dat[[gene]])
    if (sum(is.na(bn.dat)) > 0)
      next

    # construct network
    structure.1 <- hc(bn.dat)
    if (!gene %in% structure.1$arcs)
      next

    from <- c()
    to <- c()
    sub.res <- res
    sub.res[sub.res == gene] <- NA
    for (i in 1:nrow(structure.1$arcs)){
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
    training <- bn.fit(structure.2, sub.exp, ..., replace.unidentifiable=TRUE)
    models <- list.append(models, training)
    genes <- c(genes, gene)
  }
  names(models) <- genes
  return(models)
}

