setwd("C:/Users/BIGLAB_KT/Desktop/publication in work/RawData_DrCho_20221128/")

library(ggplot2)

## 1. Dr.Choi target gene validation
inputD <- read.delim("targetgene_validation.txt", header=T, sep="\t", check.names = F)
head(inputD)

inputD$status <- factor(inputD$status, levels= c("siCon","siTXNIP"))
inputD$gene <- factor(inputD$gene, levels= c("CD22","L1CAM","OTULINL","PRR5L","SDC3"))
p <- ggplot(inputD, aes(x=gene, y=value, fill=status)) + 
  geom_dotplot(binaxis='y', stackdir = "center",position = position_dodge(0.8))
# Use geom_errorbar()
p + aes(color= status)+ stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", width=0.2) +
  stat_summary(fun.y=mean, geom="point")+
  scale_y_continuous(limits = c(0,1.5), breaks = seq(0,1.5,by= 0.5))+ theme_bw()

ggsave("targetgene_validation.pdf", width = 7, height = 3, units="in")

ttest_pval <- function(inputD, gene_name){
  gene_con <- subset(inputD, gene==gene_name & status=="siCon")$value
  gene_sitxnip <- subset(inputD, gene==gene_name & status=="siTXNIP")$value
  gene_ttest <- t.test(gene_con, gene_sitxnip, paired=T)
  gene_pval <- gene_ttest$p.value
  return(gene_pval)
}

cd22_pval <- ttest_pval(inputD, "CD22")
l1cam_pval <- ttest_pval(inputD, "L1CAM")
otulinl_pval <- ttest_pval(inputD, "OTULINL")
prr5l_pval <- ttest_pval(inputD, "PRR5L")
sdc3_pval <- ttest_pval(inputD, "SDC3")

cd22_pval
l1cam_pval
prr5l_pval
otulinl_pval
sdc3_pval
#> cd22_pval [1] 0.02216107
#> l1cam_pval [1] 0.0327224
#> prr5l_pval [1] 0.6586107
#> otulinl_pval[1] 0.2951571
#> sdc3_pval [1] 0.3786222


## 2. Dr.Choi ChipAssay

inputD <- read.delim("L1CAM.chipAssay.txt", header=T, sep="\t", check.names = F)
head(inputD)

inputD$condition <- factor(inputD$condition, levels= c("siNeg","siTXNIP"))
inputD$primer <- factor(inputD$primer, levels= c("primer_3","primer_1","primer_4"))
inputD$experiment <- factor(inputD$experiment, levels= c("hdac2","H3ac"))
inputD$antibody<- factor(inputD$antibody, levels= c("target","IgG"))
p <- ggplot(inputD, aes(x=antibody, y=value, fill=condition)) + 
  geom_dotplot(binaxis='y', stackdir = "center",position = position_dodge(0.8))
# Use geom_errorbar()
p + aes(color= condition)+ stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", width=0.2) +
  stat_summary(fun.y=mean, geom="point")+ theme_bw()+
  facet_wrap(experiment ~ primer, scales = "free_y")

ggsave("L1CAM.chipAssay.pdf", width = 7, height = 4.5, units="in")
