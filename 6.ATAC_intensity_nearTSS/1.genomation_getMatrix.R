## Written by KyungTae Lee
## Last update : 2022-02-01
## R script : "1.genomation_getMatrix.R"
## Script to obtain matrix of ATAC-seq intensities (raw read counts) near TSS of expressed genes

#BiocManager::install("genomation")
library(genomation)

## test code block

#gtf.file = "Z:/reference/human/Grch38/gencode/gencode.v29.annotation.gtf"
#gtf = gffToGRanges(gtf.file)

#tss.file = "Z:/dataset/alphaArrestin_splicing/result/2021_06_03_diffBind/intersection/gene_level/gencode.v29.grch38.geneLevel.+-5kb_tss.chr18.bed"
#tss = readGeneric(tss.file, header = F, meta.col = list(geneId = 4))

#test_bam1 <- "Z:/dataset/alphaArrestin_splicing/result/2021_07_05_enrichmentHeatmap/test/bams/N48_1.trim.srt.nodup.no_chrM_MT.chr18.bam"
#test_bam2 <- "Z:/dataset/alphaArrestin_splicing/result/2021_07_05_enrichmentHeatmap/test/bams/T3-48_1.trim.srt.nodup.no_chrM_MT.chr18.bam"

#sm <- ScoreMatrix(target=test_bam1, window=tss)
#bam.files <- c(test_bam1, test_bam2)
#sm = ScoreMatrixList(targets = bam.files, windows = tss)
#sm2= ScoreMatrixBin(target = test_bam1, windows = tss, bin.num = 30)
#sm3= ScoreMatrixBin(target = test_bam2, windows = tss, bin.num = 30)
#sms <- c(sm2,sm3)
#merged <- ScoreMatrixList(targets = sms)

#multiHeatMatrix(merged, xcoords = c(-500, 500))

#heatMatrix(merged, xcoords = c(-1000, 1000))
#write.table(sm@.Data, file="C:/Users/BIGLAB_KT/Desktop/analysis/alpha_arrestin/genomation/test/test_data.txt", row.names = T, col.names = NA, sep="\t")

#pdf("/home/kyungtae/dataset/alphaArrestin_splicing/result/2021_07_05_enrichmentHeatmap/test/test.pdf")
#heatMatrix(sm, xcoords = c(-1000, 1000))
#dev.off()

## test code block ends

writeMatrix <- function(tss, inputbam, output_table){
  sm= ScoreMatrix(target = inputbam, windows = tss)
  write.table(sm@.Data, file=output_table, row.names = T, col.names = NA, sep="\t")
}

args <- commandArgs(trailingOnly = TRUE)
if (args[1] == '-h' | args[1] == '--help'){
  print ('Rscript genomation_getMatrix.R <bamdir> <tss_bed> <outputdir>')
} else {
  inputdir <- args[1]
  tss_bed <- args[2]
  outputdir <- args[3]
  inputfiles <- list.files(inputdir)
  
  tss = readGeneric(tss_bed, header = F, meta.col = list(geneId = 4))
  for (i_file in inputfiles){
    if (grepl(".bam$", i_file)){ 
      inputbam <- paste(inputdir, i_file, sep = "")
      prefix= strsplit(i_file, split=".bam")[[1]]
      outputname= paste(outputdir, prefix, ".scoreMatrix.txt", sep = "")
      writeMatrix(tss, inputbam, outputname)
    } else{}
  }
}

