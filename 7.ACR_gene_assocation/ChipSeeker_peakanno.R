## Written by KyungTae Lee
## Last update : 2021-12-02
## R script : "ChipSeeker_peakanno.R"
## Script to annotate ACR (ATAC-seq consensus peak) of genomic posiitions (promoter, gene body, downtream, intergenic, etc) by ChipSeeker R package


#BiocManager::install("ChIPseeker")
setwd("Z:/dataset/alphaArrestin_splicing/result/2021_06_03_diffBind/peakAnnotation/ChipSeeker/up5kb_down500bp/")
library(ChIPseeker)
library(GenomicFeatures)
gtf <- "Z:/reference/human/Grch38/gencode/gencode.v29.annotation.gtf"
peak_file <- "diffbind_allpeakset.bed"

txdb <- makeTxDbFromGFF(gtf)
peak <- readPeakFile(peak_file)

#options(ChIPseeker.downstreamDistance = 1000)
## downstream -> 3kb by default
peakAnno.edb <- annotatePeak(peak, tssRegion=c(-5000, 500),
                             TxDb=txdb)



peakanno_table <- as.data.frame(as.GRanges(peakAnno.edb))

prefix<- strsplit(peak_file, split=".bed")[[1]]

write.table(peakanno_table, file= paste(prefix,".ChIPseeker.tss_up5kb_down500bp.txt", sep = ""),  sep="\t", quote=F, row.names = F)

pdf(paste(prefix,".ChIPseeker.tss_up5kb_down500bp.barChart.pdf", sep = ""))
plotAnnoBar(peakAnno.edb)
dev.off()

pdf(paste(prefix,".ChIPseeker.tss_up5kb_down500bp.pieChart.pdf", sep = ""))
plotAnnoPie(peakAnno.edb)
dev.off()

pdf(paste(prefix,".ChIPseeker.tss_up5kb_down500bp.upsetplot.pdf", sep = ""))
upsetplot(peakAnno.edb)
dev.off()

