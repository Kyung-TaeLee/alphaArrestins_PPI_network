setwd("C:/Users/BIGLAB_KT/Desktop/publication in work/RawData_DrCho_20221128")

inputD <- read.delim("L1CAM.chipAssay_ttest.txt", header=T, sep="\t", check.names = F)
head(inputD)
primers <- unique(inputD$primer)
conditions <- unique(inputD$condition)
experiments <- unique(inputD$experiment)

outputdf <- data.frame("primer"=c(),"experiment"=c(),"pvalue"=c())
for (i_primer in primers) {
  for (i_experiment in experiments){
    sicon_df <- subset(inputD, primer==i_primer & experiment == i_experiment & condition=="siNeg")
    sicon_values <- sicon_df$value
    sitxnip_df <- subset(inputD, primer==i_primer & experiment == i_experiment & condition=="siTXNIP")
    sitxnip_values <- sitxnip_df$value
    res <- t.test(sicon_values, sitxnip_values, var.equal=TRUE, paired = T)
    ttest_pval <- res$p.value
    i_outputdf <- data.frame("primer"=i_primer,"experiment"=i_experiment,"pvalue"=ttest_pval)
    outputdf <- rbind(outputdf, i_outputdf)
  }
}

outputdf
write.table(outputdf, file="L1CAM.chipAssay.paired_twosided_tTest.txt", row.names = F, quote=F)


