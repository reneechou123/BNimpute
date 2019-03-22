BNimpute
<img src="../assets/logo_2.png" height="180" align="right" />
=============
Workflow:

1. Split the expression datasets into different tissue types <br>
   <b>Input:</b> 1. metadata 2. expression data <br>
   <b>Output:</b> tissue-specific epxression datasets <br>

2. For each tissue-specific expression dataset, detect outliers using WGCNA,
   remove batch effect, and then correct for PC variances. <br> 
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
neuroblast note:
#### # samples: 34
#### # genes: 17471
#### # studies: 4
``` r
data <- read.table('mushroom body_exp.tsv', header=T, row.names=1) 
datExpr0 <- as.data.frame(t(data)) 
names(datExpr0) <- rownames(data) 
rownames(datExpr0) <- names(data) 
```
#### # data cleaning
#### # removing genes not expressed in the tissue (minFraction = 0.5)
#### # excluding 1213 genes due to too many missing samples or zero variance
#### # genes: 16258
#### # cluster samples and remove outliers
#### # samples: 31
#### # plot pca before batch effect removal
``` r
data = read.table('mushroom body_cleaned.tsv', header=TRUE, row.names=1)
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
       xlim(-190, 180) +
       ylim(-190, 180)
g
```
#### # remove batch effect and plot PCA again
``` r
batch <- studies
mod.data <- data
modcombat <- model.matrix(~1, data=batch)
combat.data <- ComBat(dat=mod.data, batch=batch, mod=modcombat, par.prior=TRUE, BPPARAM=SerialParam())
combat.data.pca <- parallelPCA(combat.data, value='pca', BPPARAM=SerialParam())
var1 <- round(attr(combat.data.pca,"percentVar")[1],2)*100
var2 <- round(attr(combat.data.pca,"percentVar")[2],2)*100
combat.data.pca <- as.data.frame(combat.data.pca)
g2 <- ggplot(combat.data.pca, aes(x=-PC1, y=PC2, color=studies)) +
       geom_point() +
       scale_color_discrete(name='') +
       theme_bw() +
       theme(legend.position='bottom', legend.direction='horizontal') +
       xlab(paste0('PC1 (', var1, '%)')) +
       ylab(paste0('PC2 (', var2, '%)')) +
       xlim(-190, 180) +
       ylim(-190, 180)
g2
```
#### # correct for PC variances
``` r
n.pc <- num.sv(data, mod=modcombat, method='be', seed=123)
res <- sva_network(combat.data, n.pc)
res.pca <- parallelPCA(res, value='pca', BPPARAM=SerialParam())
var1 <- round(attr(res.pca,"percentVar")[1],2)*100
var2 <- round(attr(res.pca,"percentVar")[2],2)*100
res.pca <- as.data.frame(res.pca)
g3 <- ggplot(res.pca, aes(x=-PC1, y=PC2, color=studies)) +
       geom_point() +
       scale_color_discrete(name='') +
       theme_bw() +
       theme(legend.position='bottom', legend.direction='horizontal') +
       xlab(paste0('PC1 (', var1, '%)')) +
       ylab(paste0('PC2 (', var2, '%)')) +
       xlim(-150, 100) +
       ylim(-150, 100)
g3

data = read.table('mushroom body_cleaned.tsv', header=TRUE, row.names=1)
datExpr <- combat.data
datExpr <- as.data.frame(t(datExpr))
colnames(datExpr) <- colnames(data)
rownames(datExpr) <- rownames(data)
write.table(datExpr, 'mushroom body_cleaned_2.tsv', col.names=T, row.names=T)

res.normalized <- normalize.quantiles(res)
data = read.table('eye-antennal disc_cleaned.tsv', header=TRUE, row.names=1)
datExpr <- res.normalized
datExpr <- as.data.frame(t(datExpr))
colnames(datExpr) <- colnames(data)
rownames(datExpr) <- rownames(data)
write.table(datExpr, 'eye-antennal disc_cleaned_2_2.tsv', col.names=T, row.names=T)
```
#### # WGCNA analysis
#### # select power: 5

#### # GENIE3 analysis

#### # Cross Validation








