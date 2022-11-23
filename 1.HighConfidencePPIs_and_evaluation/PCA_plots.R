## Written by KyungTae Lee
## Last update : 2022-02-01
## R script : "PCA_plots.R"
## Script to draw PCA plots based on log2 spectral count of filtered PPIs

library(factoextra)
library(pca3d)
library(rgl)
library(plyr)

setwd("Z:/dataset/UW_protein/dataSummary_result/2017_03_08_Drosophila_Human/results_2022_01_19/pca/final_filtered/")

inputfile <- "human_spectralcount.filtered_preys.txt"
## File containing spectral counts of filtered PPI 
inputD <- read.delim(inputfile, header=T, sep="\t", check.names = F, row.names = 1)

output_prefix <- strsplit(inputfile, split = ".txt")[[1]]

pca_input <-t(log2(inputD+1))
#pca_input<- t(inputD)

#row.names(pca_input)
#colors<- mapvalues(row.names(pca_input),
#                   from = unique(row.names(pca_input)),
#                   to = c('#a6cee3', '#a6cee3','#1f78b4','#1f78b4','#1f78b4','#1f78b4','#1f78b4','#1f78b4','#1f78b4','#1f78b4','#1f78b4','#1f78b4','#1f78b4', '#33a02c','#33a02c','#fb9a99','#fb9a99', '#e31a1c','#e31a1c', '#fdbf6f','#fdbf6f',"#ff7f00","#ff7f00","#cab2d6","#cab2d6","#6a3d9a","#6a3d9a","#ffff99","#ffff99","#b15928","#b15928","grey","grey","#b2df8a","#b2df8a"))

#group <- factor(c("CG2993","CG2993","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","CG10086","CG10086","CG18748","CG18748","CG1105","CG1105","CG4674","CG4674","CG2641","CG2641","CG14696","CG14696","CG18745","CG18745","CG18744","CG18744","CG18747","CG18747","CG7047","CG7047","CG18746","CG18746"))

res.pca <- prcomp(pca_input, scale=F)
#pca3d(res.pca, group = group)
#p <- predict(res.pca)
#plot3d(p[,1:3],size = 8, type = 'p', col=colors) +
#  legend3d("topright",legend = row.names(pca_input), pch = 4, cex=0.5)

#plot3d(p[,1:3])

all_nobr_scree <- paste0(output_prefix, ".log2spec.scree.pdf", collapse = "")

fviz_eig(res.pca)
ggsave(all_nobr_scree)

ind <- get_pca_ind(res.pca )
coords <- ind$coord[,1:2]
all_nobr_pca <- paste0(output_prefix, ".log2spec.PCA.pdf", collapse = "")
## PCA plot of log2 spectral counts of filtered PPI

options(ggrepel.max.overlaps = Inf)
fviz_pca_ind(res.pca, repel = TRUE)

ggsave(all_nobr_pca, width = 7,height = 7)
