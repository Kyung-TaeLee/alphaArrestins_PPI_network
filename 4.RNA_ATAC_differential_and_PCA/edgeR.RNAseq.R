## Written by KyungTae Lee
## Last update : 2021-03-01
## R script : "edgeR.RNAseq.R"
## Script to perform analysis of differentially expressed genes in RNA-seq by edgeR and plot PCA

getGeneSymbol <- function(gene_idC){
	symbolC <- c()
	for (i_id in gene_idC){
#		ensembl_id <- strsplit(i_id, "[.]")[[1]][[1]]

		symbolC <- c(symbolC, ensembl_id)
	}
	return (symbolC)
}

writeTable<- function(colname, outputname, indataframe, genesymbol) {
	## New columna name should be given
  ## genesymbol is data.frame containing gene ID as row name and gene symbol in the first column
	
  gene_idC <- row.names(indataframe)
	#gene_symbolC <- getGeneSymbol(gene_idC)
	#gene_symbolC <- indataframe$Symbol
  gene_symbolC <- genesymbol[gene_idC,]
    writedf <- data.frame(gene_idC, gene_symbolC, indataframe)
    colnames(writedf) <- colname

    write.table(writedf, file = outputname, quote = FALSE, sep = '\t', col.names = TRUE, row.names = FALSE)
}


edger <- function(readCountTable, tableInfoFile, outPrefix, gradetreat, gradebak, l2fc_cutoff, fdr_cutoff, factor_flag)   {
  ## tableInfoFile should contain sampleId Condition Factor
  library('edgeR', verbose = F)
  library(ggplot2, verbose = F)
  library(reshape2, verbose = F)
  if (!require('hash')) {
	install.packages('hash')
  }
  library('hash')
  library(pheatmap)
  library(transcripTools)
  library(factoextra)
  library(ggrastr)
  ## For debugging
  #setwd("Z:/dataset/alphaArrestin_splicing/result/2021_05_27_TXNIP_deg/edger_batch/forRasterize")
  #readCountTable <- "TXNIP_KD3_merged_GeneExp_count.txt"
  #tableInfoFile <- "metadata.txt"
  #outPrefix <- "test"
  #gradetreat <- "siTXNIP"
  #gradebak <-"siNeg"
  #l2fc_cutoff<- 2
  #fdr_cutoff <- 0.05
  #factor_flag <- "true"
      
  l2fc_cutoff <- as.numeric(l2fc_cutoff)
  fdr_cutoff <- as.numeric(fdr_cutoff)
##read count matrix table and information table##
  readcount <- read.delim(readCountTable, header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)
  ## First column gene ID (unique), Second column gene Symbol (some redundancies may exist)

  genesymbol <- readcount[1]	  
	readcount <- readcount[c(-1)]#remove gene symbol column
	## convert to integer matrix
  readcount <- apply(readcount, c(1,2), function(x) {(as.integer(x))})
  
  ## Generation of metadata
  
  metadata  <- read.table(tableInfoFile, header = TRUE, sep = '\t', row.names = 1, check.names = FALSE)
  samples <- colnames(readcount)
  metadata <- metadata[samples,]
  group <- factor(metadata$Condition)
  group <- relevel(group, ref=gradebak)
  if (factor_flag=="true"){ 
    ## if confounding variable such as batch exist
    conf_var<- factor(metadata$Factor)
  } else{}
  
  ##DEG run##
  print ('Running edgeR')
  dge <- DGEList(counts = readcount, group= group)
  
  ## filtering out low expr
  keep <- filterByExpr(dge)
  dge <- dge[keep,,keep.lib.sizes=FALSE]
  
  ## Scaling factor for TMM
  dge <- calcNormFactors(dge) 
  
  ## MDS plot befor analysis
  mds <- paste(outPrefix,".MDS_afterTMM.pdf", sep = "")
  pdf(file=mds)
  plotMDS(dge)
  dev.off()
  
  if (factor_flag=="true"){
    design <- model.matrix(~conf_var+group)
  }else{
    design <- model.matrix(~group)
  }
  
  ## Estimating the dispersion
  dge <- estimateDisp(dge, design, robust=TRUE)
  
  ## differential expression
  fit <- glmQLFit(dge, design, robust=TRUE)
  qlf <- glmQLFTest(fit)
  
  ## Writing result table -> all
  result_n<- length(row.names(qlf$table))
  all_result <- data.frame(topTags(qlf, n=result_n))
  cpm <- cpm(dge)
  log2cpm <- log2(cpm+0.1)
  if (factor_flag=="true"){
    factor_info <- metadata[colnames(cpm),]$Factor
    condition_info <- factor(metadata[colnames(cpm),]$Condition)
	design_matrix <- model.matrix(~condition_info)
    br_log2cpm <- removeBatchEffect(log2cpm, batch=factor_info, design= design_matrix)
    br_colnames <- paste("br_log2_",colnames(br_log2cpm),sep = "")
    colnames(br_log2cpm)<- br_colnames
    cpm <- cbind(cpm, br_log2cpm[row.names(cpm),])
    
  } else {}
  all_result <- cbind(all_result, cpm[row.names(all_result),])
  newcols <- c("geneId","geneSymbol",  colnames(all_result))
  
  allresult_output <- paste(outPrefix,".allResult.txt", sep="")
  writeTable(newcols, allresult_output, all_result, genesymbol)
  
  ## Writing result table -> FDR and log2foldchange cutoff##
  head(all_result)
  sig_down <- data.frame(subset(all_result, FDR < fdr_cutoff & logFC <= -l2fc_cutoff))
  sig_up <- data.frame(subset(all_result, FDR < fdr_cutoff & logFC >= l2fc_cutoff))
  
  sig_downF <- paste0(outPrefix, ".sigDown.FDR_", fdr_cutoff, ".log2fold_", l2fc_cutoff,".txt", collapse = '')
  sig_upF <- paste0(outPrefix, ".sigUp.FDR_", fdr_cutoff, ".log2fold_", l2fc_cutoff,".txt", collapse = '')
  
  writeTable(newcols, sig_downF, sig_down, genesymbol)
  writeTable(newcols, sig_upF, sig_up, genesymbol)
  
  ## MD plot
  dt <- decideTests(qlf, adjust.method="BH", p.value = fdr_cutoff, lfc=l2fc_cutoff)
  md <- paste0(outPrefix, ".MDplot.sigDEG_log2fc_", l2fc_cutoff, "_FDR_",fdr_cutoff, ".pdf", collapse = '')
  pdf(file=md)
  #plotMD(qlf, adjust.method="BH", p.value=0.05)
  plotMD(qlf, status=dt)
  abline(h=c(-l2fc_cutoff,l2fc_cutoff), col="green")
  dev.off()

  ##volcano plot
  all_result <- na.omit(all_result)
  all_result$diffExp <- "NO"
  ## up-regulated
  all_result$diffExp[all_result$logFC > 1.0 & all_result$FDR < 0.05] <- "UP"
  ## down-regulated
  all_result$diffExp[all_result$logFC < -1.0 & all_result$FDR < 0.05] <- "DOWN"
  mycolors <- c("#2B83BA", "#D7191C", "black")
  names(mycolors) <- c("DOWN", "UP", "NO")
  
  #sigDeg <- as.factor(abs(all_result$logFC) >= as.numeric(l2fc_cutoff) & all_result$FDR <= as.numeric(fdr_cutoff))
  lgfcMax <- max(c(abs(min(all_result$logFC)-0.5), max(all_result$logFC)+0.5))
  plt <- ggplot(all_result, aes(logFC, -log10(FDR), color=diffExp)) + 
    geom_point(alpha=1, size=1) + 
    scale_color_manual(values= mycolors)+
    labs(x="log fold change", y="-log10 FDR", title=paste(c(gradetreat, 'vs', gradebak), collapse = ' ')) + 
      theme_bw() + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5),
            text = element_text(size = 16, colour = 'black'), legend.position="none")
  rasterize(plt, layerts="point", dpi=300)
  ## rasterizing
  outputVolPdf <- paste0(outPrefix, ".volcanoplot.sigDEG_log2fc_", l2fc_cutoff, "_FDR_",fdr_cutoff, ".pdf", collapse = '')
  ggsave(outputVolPdf, units = 'cm', height = 10, width = 10)
  
  ## log2CPM heatmap of degs
  sig_degs<-row.names(data.frame(subset(all_result, FDR < fdr_cutoff & abs(logFC) >= abs(l2fc_cutoff))))
  ## log2 TMMs
  
  samplenum<- length(colnames(log2cpm))
  
  
  if (factor_flag=="true"){
    br_start <- samplenum+1
    br_end <- samplenum*2
    deg_log2cpm <- cpm[sig_degs,br_start:br_end]   
    #removeBatchEffect(deg_log2cpm, batch= factor_info)
    
    
  } else{
    deg_log2cpm <- as.matrix(log2cpm[sig_degs,1:samplenum])
    }
  
  heatmap_deg <- paste0(outPrefix, ".heatmap.sigDEG_log2fc_", l2fc_cutoff, "_FDR_",fdr_cutoff, ".pdf", collapse = '')
  pdf(heatmap_deg, width=5, height=8)
  ph<- pheatmap(deg_log2cpm ,cluster_rows = T,  cluster_cols = T,
                clustering_distance_rows = "correlation",
                clustering_distance_cols = "correlation", 
                clustering_method = "ward.D",scale = "row",
                fontsize = 7, border_color = NA)
  dev.off()
    
  ## PCA analysis based on gene xpression
  ## non batch corrected CPM
  
  #log2cpm_all <- log2(cpm[,1:samplenum]+0.1)
  
  ## Selecting top 2000 variable features
  cpm_all<- log2cpm[,1:samplenum]
  top2000_var<- mostVar(data= cpm_all, n=2000)
  #top2000_var <- log2(top2000_var+0.1)
  t_top2000_var<- t(top2000_var)
  #pca_cpm_all <- t(cpm_all)
  res.pca <- prcomp(t_top2000_var, scale=F)
  all_nobr_scree <- paste0(outPrefix, ".log2CPM_all.scree.pdf", collapse = "")
  #pdf(file=all_nobr_scree)
  fviz_eig(res.pca)
  ggsave(all_nobr_scree)
  #dev.off()
  ind <- get_pca_ind(res.pca )
  coords <- ind$coord[,1:2]
  max_coord <- max(coords)
  min_coord <- min(coords)
  all_nobr_pca <- paste0(outPrefix, ".log2CPM_all.PCA.pdf", collapse = "")
  #pdf(file=all_nobr_pca)
  fviz_pca_ind(res.pca, repel = TRUE)+ xlim(min_coord, max_coord)+ ylim(min_coord, max_coord)     # Avoid text overlapping)
  #dev.off()
  ggsave(all_nobr_pca, width = 7, height = 7)
  
  ## batch corrected counts of peaksets
  #br_top2000_var <- removeBatchEffect(top2000_var, batch= factor_info)
  if (factor_flag=="true"){
    #br_top2000_var <- removeBatchEffect(top2000_var, batch= factor_info)
    #t_br_top2000_var<- t(br_top2000_var)
    #br_top2000_var<- log2(cpm[row.names(top2000_var),br_start:br_end]+0.1)
    #log2cpm_br <- log2(cpm[,br_start:br_end]+0.1)
    cpm_br_all<- cpm[,br_start:br_end]
    br_top2000_var<- mostVar(data= cpm_br_all, n=2000)
    t_br_top2000_var<- t(br_top2000_var)
    #pca_cpm_br_all <- t(cpm_br_all)
    res.pca <- prcomp(t_br_top2000_var, scale=F)
    all_br_scree <- paste0(outPrefix, ".log2CPM_all.batchRemoved.scree.pdf", collapse = "")
    #pdf(file=all_br_scree)
    fviz_eig(res.pca)
    ggsave(all_br_scree)
    #dev.off()
    ind <- get_pca_ind(res.pca )
    coords <- ind$coord[,1:2]
    max_coord <- max(coords)
    min_coord <- min(coords)
    all_br_pca <- paste0(outPrefix, ".log2CPM_all.batchRemoved.PCA.pdf", collapse = "")
    #pdf(file=all_br_pca)
    fviz_pca_ind(res.pca, repel = TRUE)+ xlim(min_coord, max_coord)+ ylim(min_coord, max_coord)
    ggsave(all_br_pca, width=7, height=7)
  #dev.off()
  } else{}
}

args <- commandArgs(trailingOnly = TRUE)
if (args[1] == '-h' | args[1] == '--help'){
    print ('Rscript edgeR.R <read count> <metadata> <output_prefix> <group of interest> <background group> <log2fold_cutoff> <FDR_cutoff> <factor_flag>')
} else {
    edger(args[1], args[2], args[3], args[4], args[5], args[6], args[7], args[8])
}

