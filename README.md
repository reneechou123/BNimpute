BNimpute
<img src="../assets/logo_2.png" height="180" align="right" />
=============
BNimpute is an R package which constructs Bayesian models with external reference RNA-seq data retrieved from the NCBI Sequence Read Archive (SRA) to impute missing values in the experimental RNA-seq count matrix. The Bayesian models are constructed using gene co-expression network (GCN) analyzed by WGCNA. The package also contains pre-trained GCN based-Bayesian imputation models for users to download and estimate the missing vlaues in <i>Drosophila</i> tissue-specific RNA-seq datasets. 


## Workflow for model construction

#### 1. Clean the tissue-specific gene expression data by following the instructions in https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-02-networkConstr-auto.pdf

#### 2. Remove batch effects for each tissue-specific gene expression data, as specified in `remove_batch_effects.R`

#### 3. Detect co-expression modules using WGCNA and retrieve the eigengene expression dataset as well as genes in each module.




```r
adjacencyA1 = adjacency(data,power=7,type="signed hybrid")
```

``` r
data <- read.table('mushroom body_exp.tsv', header=T, row.names=1) 
datExpr0 <- as.data.frame(t(data)) 
names(datExpr0) <- rownames(data) 
rownames(datExpr0) <- names(data) 
```
#### # WGCNA analysis
``` r
sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = 'signed hybrid', verbose = 5)
```
``` r
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
```
#### # select power: 6
``` r
net = blockwiseModules(datExpr, power = 7, TOMType = "signed", minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "mushboombodyTOM", networkType = 'signed hybrid', verbose = 3)
```
``` r
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],"Module colors",dendroLabels = FALSE, hang = 0.03,addGuide = TRUE, guideHang = 0.05)
```
#### # GENIE3 analysis

#### # Cross Validation








