BNimpute
<img src="../assets/logo_2.png" height="180" align="right" />
=============
Workflow:

1. Split the expression datasets into different tissue types <br>
   <b>Input:</b> 1. metadata 2. expression data <br>
   <b>Output:</b> tissue-specific epxression datasets <br>

2. For each tissue-specific expression dataset, detect outliers using WGCNA,
   drop out singleton samples for batch effect removal, remove batch effect,
   and then correct for PCA variances. <br> 
   <b>Input:</b> <br>
   <b>Output:</b> <br>
