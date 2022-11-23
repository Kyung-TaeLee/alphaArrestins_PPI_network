## Written by KyungTae Lee
## Last update : 2022-10-21
## R script : "CDF_KS_test.R"
## Script to plot CDF of log2 fold changes of gene expression or chromatin accessiblity (log2 siTXNIP/siCon)


setwd("Z:/dataset/alphaArrestin_splicing/result/2021_06_03_diffBind/acs_tf_gene_interaction/geneL2fcByAcs_denstiy_meanL2fc/up5kb_down500bp/l2fc1_fdr0.05_considered")


library(ggplot2)
library(ggpubr)
library(egg)
inputfile <- "chipseeker.tss_up5kb_down500bp.peaks_in_promoter.peak_meanl2fc_density.genebased.txt"
inputD<- read.delim(inputfile, header=T, sep="\t", check.names = F)
colnames(inputD)<- c("l2fc", "fold_flag")
head(inputD)

up_group <- subset(inputD, fold_flag=="up")$l2fc
down_group <- subset(inputD, fold_flag=="down")$l2fc
none_group <- subset(inputD, fold_flag=="nochange")$l2fc

up_group_count <- length(up_group)
down_group_count <- length(down_group)
none_group_count <- length(none_group)
  
#hdac1_group_count <- length(hdac1_group)
#hdac2_group_count <- length(hdac2_group)
#none_group_count <- length(none_group)

down_test <- ks.test(down_group, none_group, alternative = "greater")
down_pval <- down_test$p.value

up_test <- ks.test(up_group, none_group, alternative = "less")
up_pval <- up_test$p.value

down_pval_anno <- paste0("down_KStest_pval (one-sided) : ",down_pval, collapse = "")
up_pval_anno <- paste0("up_KStest_pval (one-sided) : ",up_pval, collapse = "")

text_anno <- paste(down_pval_anno, up_pval_anno, sep ="\n")

none_label <- paste0("none (N=", none_group_count, ")",collapse = "")
down_label <- paste0("closed (N=", down_group_count, ")",collapse = "")
up_label <- paste0("open (N=", up_group_count, ")",collapse = "")

#hdac1_test <- ks.test(hdac1_group, none_group, alternative = "two.sided")
#hdac1_pval <- hdac1_test$p.value

#hdac2_test <- ks.test(hdac2_group, none_group, alternative = "two.sided")
#hdac2_pval <- hdac1_test$p.value

#hdac1_pval_anno <- paste0("HDAC1_KStest_pval (two-sided) : ",hdac1_pval, collapse = "")
#hdac2_pval_anno <- paste0("HDAC2_KStest_pval (two-sided) : ",hdac2_pval, collapse = "")

#text_anno <- paste(hdac1_pval_anno, hdac2_pval_anno, sep ="\n")

#none_label <- paste0("none (N=", none_group_count, ")",collapse = "")
#hdac1_label <- paste0("HDAC1 (N=", hdac1_group_count, ")",collapse = "")
#hdac2_label <- paste0("HDAC2 (N=", hdac2_group_count, ")",collapse = "")

#p<- ggplot(inputD, aes(x=gene_l2fc, color=peak_group))+ stat_ecdf()+ coord_cartesian(xlim=c(-2,2))+ annotate(geom="text", x=-0.5, y=0.9, label= text_anno, color="black",size=1.5)+
#  theme(legend.position="top", legend.title = element_text(size=5), 
#        legend.text = element_text(size=5), 
#        legend.key.height = unit(0.5,"cm"), 
#        legend.key.width = unit(0.5,"cm"))+
#  scale_color_manual(#name="tss+-500bp_peak_based(maxMean)",
                    #breaks=c("none", "up_reg", "down_reg"),
#                    breaks= c("nochange","high","low"),
#                    labels=c(none_label, up_label, down_label),
#                    values= c("black","#d7191c","#2b83ba"))

p<- ggplot(inputD, aes(x=l2fc, color=fold_flag))+ stat_ecdf()+
  coord_cartesian(xlim=c(-2.0,2.0))+ 
  annotate(geom="text", x=-0.5, y=0.9, label= text_anno, color="black",size=1.5)+
  theme(legend.position="top", legend.title = element_text(size=5), 
        legend.text = element_text(size=5), 
        legend.key.height = unit(0.5,"cm"), 
        legend.key.width = unit(0.5,"cm"))+
    scale_color_manual(#name="tss+-500bp_peak_based(maxMean)",
    breaks=c("nochange", "up", "down"),
    #breaks= c("none","HDAC1","HDAC2"),
    labels=c(none_label, up_label, down_label),
    values= c("black","#d7191c","#2b83ba"))
p

file_prefix <- strsplit(inputfile, split=".txt")[[1]]
output<- paste(file_prefix, ".pdf", sep = "")  
ggsave(filename= output, plot = egg::set_panel_size(p=p, width=unit(7, "cm"), height=unit(7, "cm")), width=9, height = 9, units = "cm")


#tss+-5kb_peak_based(maxMean)

#hdac1_l2fc <- subset(inputD, group=="HDAC1")$l2fc
#hdac2_l2fc <- subset(inputD, group=="HDAC2")$l2fc
#none_l2fc <- subset(inputD, group=="none")$l2fc

#ks.test(hdac1_l2fc, hdac2_l2fc, alternative = "two.sided")
#ks.test(hdac1_l2fc, none_l2fc, alternative = "two.sided")
#ks.test(hdac2_l2fc, none_l2fc, alternative = "two.sided")

