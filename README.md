# BNimpute

Workflow:

1. Split the expression datasets into different tissue types <br>
   Input: 1. metadata 2. expression data <br>
   Output: tissue-specific epxression datasets <br>

2. For each tissue-specific expression dataset, detect outliers using WGCNA,
   drop out singleton samples for batch effect removal, remove batch effect,
   and then correct for PCA variances. <br> 
   Input: <br>
   Output: <br>
