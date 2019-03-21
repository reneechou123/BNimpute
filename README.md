BNimpute
<img src="../assets/logo_2.png" height="180" align="right" />
=============
Workflow:

1. Split the expression datasets into different tissue types <br>
   <b>Input:</b> 1. metadata 2. expression data <br>
   <b>Output:</b> tissue-specific epxression datasets <br>

2. For each tissue-specific expression dataset, detect outliers using WGCNA,
   remove batch effect, and then correct for PCA variances. <br> 
   <b>Input:</b> <br>
   <b>Output:</b> <br>
#
#
#
#
#
#
#
=========================================<br>
ovary note:
#### # samples: 450
#### # genes: 17471
#### # studies: 68
``` r
data <- read.table('ovary_exp.tsv', row.names=1) 
datExpr0 <- as.data.frame(t(data)) 
names(datExpr0) <- rownames(data) 
rownames(datExpr0) <- names(data) 
```
#### # data cleaning
#### # removing genes not expressed in the tissue (minFraction = 0.5)
#### # excluding 715 genes due to too many missing samples or zero variance
#### # genes: 16756
#### # cluster samples and remove outliers (at height 6e+05)
#### # samples: 437
#### # plot pca before batch effect removal
``` r
data = read.table('exp/ovary_cleaned.tsv', header=TRUE, row.names=1)
studies <- as.factor(metadata[metadata[['sample_name']] %in% rownames(data), 'study'])
data <- log2(data + 1)
data <- t(data)
data.pca <- parallelPCA(data, value='pca', BPPARAM=SerialParam())
var1 <- round(attr(data.pca,"percentVar")[1],2)*100
var2 <- round(attr(data.pca,"percentVar")[2],2)*100
data.pca <- as.data.frame(data.pca)
g <- ggplot(data.pca, aes(x=PC1, y=PC2, color=studies)) +
       geom_point() +
       scale_color_discrete(name='') +
       theme_bw() +
       theme(legend.position='bottom', legend.direction='horizontal') +
       xlab(paste0('PC1 (', var1, '%)')) +
       ylab(paste0('PC2 (', var2, '%)')) +
       xlim(-160,400) +
       ylim(-160,400)
g
```
#### # remove batch effect and plot PCA again
```r
batch <- studies
mod.data <- data
modcombat <- model.matrix(~1, data=batch)
combat.data <- ComBat(dat=mod.data, batch=batch, mod=modcombat, par.prior=TRUE, BPPARAM=SerialParam())
combat.data.pca <- parallelPCA(combat.data, value='pca', BPPARAM=SerialParam())
var1 <- round(attr(combat.data.pca,"percentVar")[1],2)*100
var2 <- round(attr(combat.data.pca,"percentVar")[2],2)*100
combat.data.pca <- parallelPCA(combat.data, value='pca')
g2 <- ggplot(combat.data.pca, aes(x=PC1, y=PC2, color=studies)) +
       geom_point() +
       scale_color_discrete(name='') +
       theme_bw() +
       theme(legend.position='bottom', legend.direction='horizontal') +
       xlab(paste0('PC1 (', var1, '%)')) +
       ylab(paste0('PC2 (', var2, '%)')) +
       xlim(-160,400) +
       ylim(-160,400)
g2
```
#### # correct for PCA variances

Cross Validation








