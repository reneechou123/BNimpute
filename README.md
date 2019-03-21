BNimpute
<img src="../assets/logo_2.png" height="180" align="right" />
=============
Workflow:

1. Split the expression datasets into different tissue types <br>
   <b>Input:</b> 1. metadata 2. expression data <br>
   <b>Output:</b> tissue-specific epxression datasets <br>

2. For each tissue-specific expression dataset, detect outliers using WGCNA,
   drop out singleton samples for batch effect removal, remove batch effect,
   and then correct for PCA variances. <br> datExpr0 = as.data.frame(t(data))
   <b>Input:</b> <br>
   <b>Output:</b> <br>
   
=========================================
ovary note:
#### samples: 450
#### genes: 17471
#### studies: 68
data <- read.table('ovary_exp.tsv', row.names=1)
datExpr0 <- as.data.frame(t(data))
names(datExpr0) <- rownames(data)
rownames(datExpr0) <- names(data)
#### data cleaning
#### removing genes not expressed in the tissue (minFraction = 0.5)
#### excluding 715 genes due to too many missing samples or zero variance
#### genes: 16756


Cross Validation
