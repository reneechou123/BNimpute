BNimpute
<img src="../assets/logo_2.png" height="180" align="right" />
=============
BNimpute is an R package which constructs Bayesian models with external reference RNA-seq data retrieved from the NCBI Sequence Read Archive (SRA) to impute missing values in the experimental RNA-seq count matrix. The Bayesian models are constructed using gene co-expression network (GCN) analyzed by WGCNA. The package also contains pre-trained GCN based-Bayesian imputation models for users to download and estimate the missing vlaues in <i>Drosophila</i> tissue-specific RNA-seq datasets. 


## Workflow for model construction

#### 1. Clean the tissue-specific gene expression data by following the instructions in https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf

#### 2. Remove batch effects for each tissue-specific gene expression data, as specified in `remove_batch_effects.R`

#### 3. Detect co-expression modules using WGCNA and retrieve the eigengene expression dataset as well as genes in each module, as specified in `WGCNA_analysis.R`




```r
adjacencyA1 = adjacency(data,power=7,type="signed hybrid")
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
