BNimpute
<img src="../assets/logo_2.png" height="180" align="right" />
=============
BNimpute is an R package which constructs Bayesian models with external reference RNA-seq data retrieved from the NCBI Sequence Read Archive (SRA) to impute missing values in the experimental RNA-seq count matrix. The Bayesian models are constructed using gene co-expression network (GCN) analyzed by WGCNA. The package also contains pre-trained GCN based-Bayesian imputation models for users to download and estimate the missing vlaues in <i>Drosophila</i> tissue-specific RNA-seq datasets. 


## Workflow for preparing input datasets

#### 1. Clean the tissue-specific gene expression data by following the instructions in https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf

#### 2. Remove batch effects for each tissue-specific gene expression data, as specified in `remove_batch_effects.R`

#### 3. Detect co-expression modules using WGCNA and retrieve the eigengene expression dataset as well as genes in each module, as specified in `WGCNA_analysis.R`


```r
adjacencyA1 = adjacency(data,power=7,type="signed hybrid")
```
