library(WGCNA)
library(flashClust)
library(dynamicTreeCut)


#=====================================================================================
#
#  Import data
#
#=====================================================================================

ref.exp <- read.table('mushroom body_cleaned_2.tsv', header=T, row.names=1)


#=====================================================================================
#
#  Choose the soft-threshold power
#  (Adapted from the WGCNA tutorial: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.pdf)
#
#=====================================================================================

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))

# Call the network topology analysis function
sft <- pickSoftThreshold(ref.exp, powerVector = powers, networkType = 'signed hybrid', verbose = 5)

# Plot the results
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")

# this line corresponds to using an R^2 cut-off of h
abline(h=0.95,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")


#=====================================================================================
#
#  Detect gene co-expression modules
#  (Adapted from the WGCNA tutorial: https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/JMiller/Tutorial%20document.pdf)
#
#=====================================================================================

# Build gene dendrogram
adjacency <- adjacency(ref.exp, power=8, type="signed hybrid")
diag(adjacency) <- 0
dissTOM <- 1 - TOMsimilarity(adjacency, TOMType="signed")
geneDendro <- flashClust(as.dist(dissTOM), method="average")

# Cut the dendrogram to generate co-expression modules
modules <- NULL
ds <- 0
tree <- cutreeHybrid(dendro = geneDendro, pamStage=FALSE,
		    minClusterSize = (30-3*ds), cutHeight = 0.99,
		    deepSplit = ds, distM = dissTOM)
modules <- cbind(modules, labels2colors(tree$labels))
rownames(modules) <- colnames(ref.exp)

# Plot the dendrogram with module information and calculate eigen expressions for each module
plotDendroAndColors(geneDendro, modules, "Modules", main = "", dendroLabels=FALSE)
eig.exp <- moduleEigengenes(ref.exp, colors=modules)$eigengenes


