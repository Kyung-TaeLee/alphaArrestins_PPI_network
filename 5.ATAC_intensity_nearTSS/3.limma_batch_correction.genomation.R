## Written by KyungTae Lee
## Last update : 2022-04-17
## R script : "3.limma_batch_correction.genomation.R"
## Script to correct ATAC-seq read counts from different batches of batch effect

library(limma)

args <- commandArgs(trailingOnly = TRUE)
if (args[1] == '-h' | args[1] == '--help'){
  print ('Rscript limma_batch_correction.R <geneSampleTable> <batches>')
} else {
  geneSampleTable <- args[1]
  conditions <- args[2]
  batches <- args[3]

  inputD <- read.delim(geneSampleTable, header=T, row.names = 1, check.names = F)
  condition <- factor(as.vector(strsplit(conditions, split = ",")[[1]]))
  design_matrix <- model.matrix(~condition)
  batch <- as.vector(strsplit(batches, split = ",")[[1]])
  log2input <- log2(inputD+1)
  br_matrix <- removeBatchEffect(log2input, batch= batch, design= design_matrix)
  output_prefix <- strsplit(geneSampleTable, split = ".txt")[[1]]
  outputname <- paste(output_prefix, ".log2.limma_batchCorrect.txt", sep = "")
  br_matrix <- cbind(data.frame("probe"=rownames(br_matrix)), br_matrix)
  write.table(br_matrix, file=outputname, row.names = F,sep="\t", quote = F)
}

