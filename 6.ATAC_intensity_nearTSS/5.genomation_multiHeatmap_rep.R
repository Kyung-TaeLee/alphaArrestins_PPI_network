## Written by KyungTae Lee
## Last update : 2022-04-17
## R script : "5.genomation_multiHeatmap_rep.R"
## Script to plot heatmap of batch-corrected, normalized ATAC-seq intensities near near TSS of expressed genes

library(genomation)

args <- commandArgs(trailingOnly = TRUE)
if (args[1] == '-h' | args[1] == '--help'){
  #print ('Rscript genomation_multiHeatmap_rep.R <inputdir> <order_flag> <batch_flag> <gene_l2fc_file> <bound>')
  print ('Rscript genomation_multiHeatmap_rep.R <inputdir> <order_flag> <batch_flag> <bound>')
} else {

  inputdir <- args[1]  
  order_flag <- args[2]
  batchCor_flag <- args[3]
  #gene_l2fc_file <- args[4]
  bound <- args[4]
  
  setwd(inputdir)
  if (batchCor_flag=="true"){
    prefix= ".batchCorrect"
  } else{
    prefix= ""
  }
 
  #gene_l2fc <- read.delim(gene_l2fc_file, header=T, sep="\t", check.names = F)
  
  #up_group<-which(gene_l2fc$peak_average_l2fc>=1)
  #down_group <- which(gene_l2fc$peak_average_l2fc<=-1)
  
  #target_group <- c(up_group, down_group)
  #target_genes <- gene_l2fc[target_group,]
  #write.table(target_genes$gene, file=paste0("target_gene.txt"), sep="\t", quote=F, row.names = F)
  
  bound <- as.numeric(bound)
    
  sineg1 <- read.delim(paste("siNeg_1", prefix, ".scoreMatrix.txt", sep = ""), header=T, check.names = F,sep="\t", row.names = 1)
  
  sineg1_mat<- as.matrix(sineg1)
  ## add pseudocount 1
  temp_mat <- 2^sineg1_mat
  temp_mat <- temp_mat+ 1
  sineg1_mat <- log2(temp_mat)
  
  colnames(sineg1_mat) <- -(bound+1):bound
  #sineg1_sm <- as(sineg1_mat[target_group,], "ScoreMatrix")
  sineg1_sm <- as(sineg1_mat, "ScoreMatrix")
  
  sineg2 <- read.delim(paste("siNeg_2", prefix, ".scoreMatrix.txt", sep = ""), header=T, check.names = F,sep="\t", row.names = 1)
  
  sineg2_mat<- as.matrix(sineg2)
  ## add pseudocount 1
  temp_mat <- 2^sineg2_mat
  temp_mat <- temp_mat+ 1
  sineg2_mat <- log2(temp_mat)
  
  colnames(sineg2_mat) <- -(bound+1):bound
  #sineg2_sm <- as(sineg2_mat[target_group,], "ScoreMatrix")
  sineg2_sm <- as(sineg2_mat, "ScoreMatrix")
  
  sitxnip1 <- read.delim(paste("siTXNIP_1", prefix, ".scoreMatrix.txt", sep = ""), header=T, check.names = F,sep="\t", row.names = 1)
  
  sitxnip1_mat<- as.matrix(sitxnip1)
  ## add pseudocount 1
  temp_mat <- 2^sitxnip1_mat
  temp_mat <- temp_mat+ 1
  sitxnip1_mat <- log2(temp_mat)
  
  colnames(sitxnip1_mat) <- -(bound+1):bound
  #sitxnip1_sm <- as(sitxnip1_mat[target_group,], "ScoreMatrix")
  sitxnip1_sm <- as(sitxnip1_mat, "ScoreMatrix")
  
  sitxnip2 <- read.delim(paste("siTXNIP_2", prefix, ".scoreMatrix.txt", sep = ""), header=T, check.names = F,sep="\t", row.names = 1)
  
  sitxnip2_mat<- as.matrix(sitxnip2)
  ## add pseudocount 1
  temp_mat <- 2^sitxnip2_mat
  temp_mat <- temp_mat+ 1
  sitxnip2_mat <- log2(temp_mat)
  
  colnames(sitxnip2_mat) <- -(bound+1):bound
  #sitxnip2_sm <- as(sitxnip2_mat[target_group,], "ScoreMatrix")
  sitxnip2_sm <- as(sitxnip2_mat, "ScoreMatrix")
  
  sms<- list(sineg1_sm, sineg2_sm, sitxnip1_sm, sitxnip2_sm)
  sm = ScoreMatrixList(targets = sms)
  names(sm)<- c("siNeg1","siNeg2", "siTXNIP1","siTXNIP2")
  

  #multiHeatMatrix(sm, xcoords = c(-3000, 3000), winsorize = c(0,99), common.scale = T, col = c("lightgray", "blue"), order=TRUE)
  if (order_flag=="false"){
    pdf(paste0("multiheatmap.pdf", collapse = ""))
    #multiHeatMatrix(sm, xcoords = c(-500, 500), common.scale = T, col = c("lightgray", "blue"), group=list(open= up_group, none= none_group, closed= down_group))
    #multiHeatMatrix(sm, xcoords = c(-bound, bound), common.scale = T, col = c("lightgray", "blue"), group= list(open=up_group, closed= down_group))
    multiHeatMatrix(sm, xcoords = c(-bound, bound), common.scale = T, col = c("lightgray", "blue"))
    dev.off()
    
  } else if(order_flag=="true"){
    pdf("multiheatmap.ordered.pdf")
    row_order<-multiHeatMatrix(sm, xcoords = c(-bound, bound), common.scale = T, col = c("lightgray", "blue"), order=TRUE, winsorize = c(0,99))
    dev.off()
    row_order_names <- names(row_order)
    write.table(row_order_names, file="sorted_row_order.txt", row.names = F, quote=F)
  }
  
}
