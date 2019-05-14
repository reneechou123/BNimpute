library(feather)
library(sva)
library(scran)
library(ggplot2)


#=====================================================================================
#
#  Import data
#
#=====================================================================================

data = read.table('mushroom body_cleaned.tsv', header=TRUE, row.names=1)
metadata <- as.data.frame(read_feather('rnaseq_metadata_copy_feather.tsv'))
studies <- as.factor(metadata[metadata[['sample_name']] %in% rownames(data), 'study'])
data <- log2(data + 1)
data <- t(data)


#=====================================================================================
#
#  Plot PCA before batch effect removal
#
#=====================================================================================


data.pca <- parallelPCA(data, value='pca', BPPARAM=SerialParam())
var1 <- round(attr(data.pca,"percentVar")[1],2)*100
var2 <- round(attr(data.pca,"percentVar")[2],2)*100
data.pca <- as.data.frame(data.pca)
g <- ggplot(data.pca, aes(x=PC1, y=PC2, color=studies)) +
       geom_vline(xintercept=0, linetype='dashed', size=1, alpha=0.5) +
       geom_hline(yintercept=0, linetype='dashed', size=1, alpha=0.5) +
       geom_point(size=5, alpha=0.8) +
       # scale_color_discrete() +
       scale_color_brewer(palette = "Dark2") +
       theme_bw() +
       theme(legend.position = "none",
             legend.title=element_text(size=18, face='bold'),
             legend.text=element_text(size=15),
             plot.title=element_text(size=25, face="bold", hjust=0.5),
             axis.text=element_text(size=20, face='bold'),
             axis.title=element_text(size=25,face="bold")) +
       xlab(paste0('PC1 (', var1, '%)')) +
       ylab(paste0('PC2 (', var2, '%)')) +
       ggtitle('Mushroom body') +
       xlim(min(min(data.pca$PC1), min(data.pca$PC2)), max(max(data.pca$PC1), max(data.pca$PC2))) +
       ylim(min(min(data.pca$PC1), min(data.pca$PC2)), max(max(data.pca$PC1), max(data.pca$PC2)))
g


#=====================================================================================
#
#  Remove batch effects
#
#=====================================================================================


batch <- studies
mod.data <- data
modcombat <- model.matrix(~1, data=batch)
combat.data <- ComBat(dat=mod.data, batch=batch, mod=modcombat, par.prior=TRUE, BPPARAM=SerialParam())

# plot pca after batch effect removal
combat.data.pca <- parallelPCA(combat.data, value='pca', BPPARAM=SerialParam())
var1 <- round(attr(combat.data.pca,"percentVar")[1],2)*100
var2 <- round(attr(combat.data.pca,"percentVar")[2],2)*100
combat.data.pca <- as.data.frame(combat.data.pca)
g2 <- ggplot(combat.data.pca, aes(x=-PC1, y=PC2, color=studies)) +
       geom_vline(xintercept=0, linetype='dashed', size=1, alpha=0.5) +
       geom_hline(yintercept=0, linetype='dashed', size=1, alpha=0.5) +
       geom_point(size=5, alpha=0.8) +
       scale_color_brewer(palette = "Dark2") +
       theme_bw() +
       theme(legend.position = "none",
             plot.title=element_text(size=25, face="bold", hjust=0.5),
             axis.text=element_text(size=20, face='bold'),
             axis.title=element_text(size=25,face="bold")) +
       xlab(paste0('PC1 (', var1, '%)')) +
       ylab(paste0('PC2 (', var2, '%)')) +
       xlim(min(min(data.pca$PC1),min(data.pca$PC2)) - 0.05, max(max(data.pca$PC1),max(data.pca$PC2)) + 0.05) +
       ylim(min(min(data.pca$PC1),min(data.pca$PC2)) - 0.05, max(max(data.pca$PC1),max(data.pca$PC2)) + 0.05)  
g2


#=====================================================================================
#
#  Export data
#
#=====================================================================================


data = read.table('mushroom body_cleaned.tsv', header=TRUE, row.names=1)
datExpr <- 2^(combat.data) - 1
datExpr <- as.data.frame(t(datExpr))
colnames(datExpr) <- colnames(data)
rownames(datExpr) <- rownames(data)
write.table(datExpr, 'mushroom body_cleaned_2.tsv', col.names=T, row.names=T)


