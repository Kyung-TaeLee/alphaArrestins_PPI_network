## Written by KyungTae Lee
## Last update : 2022-05-21
## R script : "orthologInteractions.clustering.R"
## Script to plot hierarchical clustered heatmap of PPIs between alpha-arrestin and ortholog protein groups between human and fly

library(pheatmap)
library(dplyr)
library(ggplot2)
library(factoextra)

setwd("Z:/dataset/UW_protein/dataSummary_result/2017_03_08_Drosophila_Human/Ortholog/Network analysis/2022_01_17/bait_clustering/")

inputD <- read.delim("ortholog_interaction.txt", header=T, sep="\t", row.names = 1, check.names = F)
anno <- read.delim("ortholog_interactiongroup_anno.txt", header=T, sep="\t", row.names = 1, check.names = F)
## outputs from getOrthoInterTable.py should be given

## remove ortholog gene group without functionanl information
#func_ortholog <- row.names(subset(anno, collapsed_goterm!="other"))
#inputD <- inputD[func_ortholog,]

input_inter <- t(log2(inputD+0.1))
anno <- data.frame("collapsed_goterm"= anno[colnames(input_inter),])
row.names(anno) <- colnames(input_inter)

#pdf("ortholog_inter.heatmap.correlation_ward.D_clustering.pdf", height = 7, width= 12)
#ph <- pheatmap(input_inter, show_rownames = T, show_colnames= T, cluster_rows = T, cluster_cols = T, clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", clustering_method = "ward.D", cutree_rows = 7, annotation_col = anno)
#dev.off()


d<- as.dist(1 - cor(t(input_inter)))
fit <- hclust(d, method="ward.D")
pdf("ortholog_inter.heatmap.correlation_ward.D_clustering.hc_tree_k7.pdf")
plot(fit, hang=-1)
rect.hclust(fit, k=7)
dev.off()

