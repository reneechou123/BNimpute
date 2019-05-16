BNimpute
<img src="../assets/logo_2.png" height="180" align="right" />
=============
BNimpute is an R package which constructs tissue-specific Bayesian gene models with external reference RNA-seq data retrieved from the NCBI Sequence Read Archive (SRA). The goal is to impute missing values in experimental RNA-seq count matrices. The Bayesian models are constructed using gene co-expression network (GCN) analyzed by WGCNA. The package also contains pre-trained GCN based-Bayesian imputation models for users to download and estimate the missing vlaues in <i>Drosophila</i> RNA-seq datasets. 


## Workflow for preparing input datasets

#### 1. Clean the tissue-specific gene expression data by following the instructions in https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/FemaleLiver-01-dataInput.pdf

#### 2. Remove batch effects for each tissue-specific gene expression data, as specified in `remove_batch_effects.R`

#### 3. Detect co-expression modules using WGCNA and retrieve modules as well as the eigengene expression dataset, as specified in `WGCNA_analysis.R`

##### parallel
##### build.model
