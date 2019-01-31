# BNimpute

Workflow:

1. Split the expression datasets into different tissue types
   Input: 1. metadata 2. expression data
   Output: tissue-specific epxression datasets

2. For each tissue-specific expression dataset, detect outliers using WGCNA,
   drop out singleton samples () 
