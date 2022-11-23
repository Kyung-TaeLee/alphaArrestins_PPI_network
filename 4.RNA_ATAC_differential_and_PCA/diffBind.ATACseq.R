## Written by KyungTae Lee
## Last update : 2022-04-05
## R script : "diffBind.ATACseq.R"
## Script to 1. obtain consensus peaks from multiple ATAC-seq samples , 2. perform differential analysis (edgeR) accessibility of chromatin regions between two conditions, 3. plot PCA

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#install.packages("Matrix")
#BiocManager::install("DiffBind")

library(DiffBind)
library(limma)
library(pheatmap)
library(transcripTools)
library(factoextra)
library(Repitools)

diffbind_analysis <-function(sample_sheet, ctrl_condition, fdr_cutoff, fold_cutoff,  output_prefix){
  ## sample_sheet should contain
  ## SampleID, Tissue, Factor (used for confounding variable), Condition (used for contrast),  bamReads, Peaks, PeakCaller

  ## for debugging
  #setwd("C:/Users/BIGLAB_KT/Desktop/analysis/alpha_arrestin/diffBind/")
  #sample_sheet <- "chr18_test.csv"
  #ctrl_condition <- "siNeg"
  #fdr_cutoff <- 0.05
  #fold_cutoff <- 2
  #output_prefix<- "test_output/test"
  
  fdr_cutoff <- as.numeric(fdr_cutoff)
  fold_cutoff <- as.numeric(fold_cutoff)  
  ## generation of dba object
  samples <- read.csv(sample_sheet)
  samples <- dba(sampleSheet =samples)
  ## Write consensus peakset
  #consensus_peak <- paste(output_prefix, ".consensusPeak.txt", sep = "")
  #consensus <- dba.peakset(samples, minOverlap = 2, bRetrieve = T, writeFile = consensus_peak)
  
  ## Counting reads 
  read_count <- dba.count(samples,minOverlap = 2)
  
  ## Normalizing data (based on sequencing depth)
  readcount_norm <- dba.normalize(read_count, method = DBA_EDGER)
  rcnorm<- dba.normalize(read_count, method = DBA_EDGER, bRetrieve = T)
  
  ## PCA based on normalized counts of consensus peaks
  rcnorm_pca <- paste(output_prefix, ".Normalized_readCount.PCA.pdf", sep="")
  pdf(file= rcnorm_pca)  
  dba.plotPCA(readcount_norm, label="ID")
  dev.off()
  
  ## contrast design
  #df_contrast <- dba.contrast(readcount_norm, minMembers = 2,
                              #reorderMeta = list(Treatment=ctrl_condition))
  df_contrast <-dba.contrast(readcount_norm, design="~ Factor + Treatment", minMembers = 2,contrast= c("Treatment","siTXNIP","siNeg"))

  ## Performing differential analysis
  test_df <- dba.analyze(df_contrast, method = DBA_EDGER)
  
  test_result <- paste0(output_prefix, ".FDR_",fdr_cutoff, ".log2fc_",fold_cutoff,".peakset.txt", collapse = "")
  ## Write down only differential peaks
  test_report <- dba.report(test_df,method=DBA_EDGER, bUsePval=F, th=fdr_cutoff,fold=fold_cutoff, bCounts = T, bNormalized = T)
  write.table(test_report, file=test_result, quote=F, sep="\t", row.names = F)
  ## Write down all peaks
  #all_result <- paste0(output_prefix, ".all.peakset.txt", collapse = "")
  all_report <- dba.report(test_df,method=DBA_EDGER, th= 1, bCounts = T, bNormalized = T)
  #write.table(all_report, file=all_result, quote=F, sep="\t", row.names = F)
  

  ## PCA based on differentially bound sites 
  diff_pca <- paste0(output_prefix, ".FDR_",fdr_cutoff, ".log2fc_",fold_cutoff,".PCA.pdf", collapse = "")
  pdf(file=diff_pca)
  dba.plotPCA(test_df, contrast=1, method=DBA_EDGER, label="Replicate", bUsePval=F, th= fdr_cutoff,fold=fold_cutoff)
  dev.off()

  ## MA plots
  diff_ma <- paste0(output_prefix, ".FDR_",fdr_cutoff, ".log2fc_",fold_cutoff,".MA.pdf", collapse = "")
  pdf(file=diff_ma)
  dba.plotMA(test_df, method=DBA_EDGER, bUsePval=F, th=fdr_cutoff,fold=fold_cutoff,bNormalized=T)
  dev.off()

  ## Volcano plots
  diff_volcano <- paste0(output_prefix, ".FDR_",fdr_cutoff, ".log2fc_",fold_cutoff,".Volcano.pdf", collapse = "")
  pdf(file=diff_volcano)
  dba.plotVolcano(test_df, method=DBA_EDGER, bUsePval=F, th=fdr_cutoff,fold=fold_cutoff)
  dev.off()

  ## Box plot of read distributions
  diff_boxplot <- paste0(output_prefix, ".FDR_",fdr_cutoff, ".log2fc_",fold_cutoff,".readcount_boxplot.pdf", collapse = "")
  pdf(file=diff_boxplot)
  pvals<- dba.plotBox(test_df, method=DBA_EDGER, bUsePval=F, th=fdr_cutoff,fold=fold_cutoff)
  dev.off()
  
  pval_output <- paste(output_prefix, ".readcount_pvals.txt", sep="")
  write.table(pvals, file=pval_output, sep="\t", quote=F)
  
  ## batch correction and write table with counts 
  all_acs_report <- annoGR2DF(all_report)
  sample_num<-length(row.names(samples$samples))
  col_num<- length(colnames(all_acs_report))
  count_range <- (col_num-sample_num+1):(col_num)
  
  norm_count_all <- all_acs_report[,count_range]
  
  ## Batch effect removed counts  
  samples
  batch_info <- data.frame("batch"= samples$samples$Factor)
  condition_info <- data.frame("condition"= samples$samples$Condition)
  row.names(batch_info) <- samples$samples$SampleID
  condition_info$condition <- factor(condition_info$condition)
  row.names(condition_info) <- samples$samples$SampleID
  count_cols <- colnames(norm_count_all)
  sorted_batchinfo <- batch_info[count_cols,]
  sorted_condition <- factor(condition_info[count_cols,])
  design_matrix <- model.matrix(~sorted_condition)
    
  log2count_all <- log2(norm_count_all+0.1)
  br_log2count_all <- removeBatchEffect(log2count_all, batch= sorted_batchinfo, design= design_matrix)
  newcolnames <- paste("br_log2_", colnames(br_log2count_all), sep="")
  colnames(br_log2count_all)<- newcolnames
  allacs_bradded <- cbind(all_acs_report,br_log2count_all)
  all_result <- paste0(output_prefix, ".all.peakset.normcount.batchCorrect_normcount.txt", collapse = "")
  write.table(allacs_bradded, file=all_result, quote=F, sep="\t", row.names = F)
  
  ## Heatmap of differential features 
  de_acs_report <- annoGR2DF(test_report)
  
  br_count_range <- count_range + sample_num
  de_count_data <- allacs_bradded[row.names(de_acs_report), br_count_range]
  #br_de_count_data <- removeBatchEffect(de_count_data, batch= sorted_batchinfo)
  #br_de_count_data[br_de_count_data<0]<-0
  #br_de_log2count <- log2(br_de_count_data +0.1)
  #heatmap_countdata <- count_data
  #heatmap_countdata[heatmap_countdata<0]<-0
  #log2count <- log(heatmap_countdata+0.1)
  
  diff_heatmap <- paste0(output_prefix, ".FDR_",fdr_cutoff, ".log2fc_",fold_cutoff,".heatmap.pdf", collapse = "")
  pdf(file=diff_heatmap)
  ph<- pheatmap(de_count_data,cluster_rows = T,  cluster_cols = T,
                clustering_distance_rows = "correlation",
                clustering_distance_cols = "correlation", 
                clustering_method = "ward.D",scale = "row",
                fontsize = 7, border_color = NA)
  dev.off()
  
  ## PCA analysis based on consensus peaksets
  ## non batch corrected counts of peaksets

  ## Selecting top 2000 variable features

  top2000_var<- mostVar(data= log2count_all, n=2000)
  t_top2000_var<- t(top2000_var)
  #pca_count_all <- t(log2(norm_count_all+0.1))
  res.pca <- prcomp(t_top2000_var, scale=F)
  all_nobr_scree <- paste0(output_prefix, ".conseusen_peakset_all.scree.pdf", collapse = "")
  #pdf(file=all_nobr_scree)
  fviz_eig(res.pca)
  ggsave(all_nobr_scree)
  #dev.off()
  ind <- get_pca_ind(res.pca )
  coords <- ind$coord[,1:2]
  max_coord <- max(coords)
  min_coord <- min(coords)
  all_nobr_pca <- paste0(output_prefix, ".conseusen_peakset_all.PCA.pdf", collapse = "")
  #pdf(file=all_nobr_pca)
  fviz_pca_ind(res.pca, repel = TRUE)+ xlim(min_coord, max_coord)+ ylim(min_coord, max_coord)     # Avoid text overlapping)
  #dev.off()
  ggsave(all_nobr_pca, width = 7,height = 7)
  
  ## batch corrected counts of peaksets
  #log2count_br <- log2(br_normcount_all+0.1)
  br_top2000_var <- mostVar(data= br_log2count_all, n=2000)
  #log2count_br <- removeBatchEffect(top2000_var, batch= sorted_batchinfo)
  t_br_top2000_var<- t(br_top2000_var)
  #pca_br_count_all <- t(log2(br_normcount_all+0.1))
  res.pca <- prcomp(t_br_top2000_var, scale=F)
  all_br_scree <- paste0(output_prefix, ".conseusen_peakset_all.batchRemoved.scree.pdf", collapse = "")
   #pdf(file=all_br_scree)
   fviz_eig(res.pca)
   #dev.off()
   ggsave(all_br_scree)
   ind <- get_pca_ind(res.pca )
   coords <- ind$coord[,1:2]
   max_coord <- max(coords)
   min_coord <- min(coords)
   all_br_pca <- paste0(output_prefix, ".conseusen_peakset_all.batchRemoved.PCA.pdf", collapse = "")
   #pdf(file=all_br_pca)
   fviz_pca_ind(res.pca, repel = TRUE)+ xlim(min_coord, max_coord)+ ylim(min_coord, max_coord)
   #dev.off()
   ggsave(all_br_pca, width = 7,height = 7)
  
}


args <- commandArgs(trailingOnly = TRUE)
if (args[1] == '-h' | args[1] == '--help'){
  print ('Rscript diffBind.R <sample_sheet> <ctrl_condition> <fdr_cutoff> <log2fc_cutoff> <output_prefix>')
} else {
  diffbind_analysis(args[1], args[2], args[3], args[4], args[5])
}
