## Written by KyungTae Lee
## Last update : 2022-02-01
## R script : "highConfidenceHeatmap.R"
## Script to plot log2 spectral count heatmap with hierarchical clustering

setwd("Z:/dataset/UW_protein/dataSummary_result/2017_03_08_Drosophila_Human/Human_protein/SAINT_analysis/SAINTexpress/output/R2_L3/bait_prey_2dTable/")
getwd()

library(ComplexHeatmap)
library(dendextend)

fly_anno <- c("Nedd4","CG42797","Smurf","atl","Su(dx)","kibra" ,"yki","HERC2")
## Known proteins to intearct with alpha-arrestins in fly in the high-confidence PPIs
human_anno <- c("ECD", "HECW2", "ITCH", "NEDD4", "NEDD4L", "TSG101","WWP1","WWP2","ZMYM2")
## Known proteins to intearct with alpha-arrestins in human in the high-confidence PPIs

input_file <- "human.filtered_preys.spectralcount_heatmap_input.txt"
## specral count formatted as 2D matrix
prefix <- strsplit(input_file, split = ".txt")[[1]]

indata = read.table(input_file ,sep='\t',header=TRUE,row.names=1, check.names = F)
indata<- round(log2(indata+1),digits=3)
## get colors for the heatmap
max_spec<- max(indata)
min_spec<- min(indata)
interval<- (max_spec - min_spec)/20
colors <- c(seq(min_spec,max_spec, by= interval))
colorLen<-length(colors)
hmcol <- colorpanel(colorLen-1, '#FBF8D7','#35ABC3' ,'#1D2141')

output_pdf <- paste(prefix,".pdf", sep="")
pdf(output_pdf)
hm<- Heatmap(indata, name = "mat", col = hmcol, clustering_distance_rows = "pearson", clustering_distance_columns = "pearson",clustering_method_rows = "ward.D", clustering_method_columns = "ward.D", column_split = 6) 
print(hm)
dev.off()
col_order <- column_order(hm)
length(col_order)
for (i_index in 1:length(col_order)){
  indexes <- col_order[[i_index]]
  gene_symbol <- colnames(indata)[indexes]
  output_name <- paste0(prefix, ".cluster_", i_index, ".txt", collapse = "")
  write.table(gene_symbol, file=output_name, quote = F, row.names = F)
}

sessionInfo()


