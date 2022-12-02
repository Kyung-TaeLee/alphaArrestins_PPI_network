setwd("C:/Users/BIGLAB_KT/Desktop/publication in work/RawDataNgoc_20221125")

library(ggplot2)


## 1. Ngoc western blot of siNegative and siTXNIP
inputD <- read.delim("westernblot.txt", header=T, sep="\t", check.names = F)

p <- ggplot(inputD, aes(x=status, y=value)) + 
  geom_dotplot(binaxis='y', stackdir='center')

# Use geom_errorbar()
p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="red", width=0.2) +
  stat_summary(fun.y=mean, geom="point", color="red")+
  scale_y_continuous(limits = c(0,1.3), breaks = seq(0,1.3,by= 0.3))+ theme_bw()
  
ggsave("westernblot.pdf")
# paired student's t-test

sicon_values <- subset(inputD, status=="siCon")$value
sitxnip_values <- subset(inputD, status=="siTXNIP")$value
ttest <- t.test(sicon_values, sitxnip_values, alternative="two.sided", paired=T)
ttest
# t = 6.3127, df = 2.4321, p-value = 0.0146

## 2. Ngoc qPCR data of siNegative and siTXNIP
inputD <- read.delim("qpcr_tubulin.txt", header=T, sep="\t", check.names = F)

p <- ggplot(inputD, aes(x=status, y=value)) + 
  geom_dotplot(binaxis='y', stackdir='center')

# Use geom_errorbar()
p + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), 
                 geom="errorbar", color="red", width=0.2) +
  stat_summary(fun.y=mean, geom="point", color="red")+
  scale_y_continuous(limits = c(0,1.2), breaks = seq(0,1.2,by= 0.3))+ theme_bw()
ggsave("qpcr_tubulin.pdf")

# paired student's t-test

sicon_values <- subset(inputD, status=="siCon")$value
sitxnip_values <- subset(inputD, status=="siTXNIP")$value
ttest <- t.test(sicon_values, sitxnip_values, alternative="two.sided", paired=T)
ttest
# t = 15.793, df = 3.3268, p-value = 0.0003411

## 3. Ngoc fractionation western blot

inputD <- read.delim("fractionation_westernblot.txt", header=T, sep="\t", check.names = F)
head(inputD)
inputD$fraction<- factor(inputD$fraction, levels= c("C","N"))
p <- ggplot(inputD, aes(x=fraction, y=value, fill=status)) + 
  geom_dotplot(binaxis='y', stackdir = "center",position = position_dodge(0.8))
# Use geom_errorbar()
p + aes(color= status)+ stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="errorbar", width=0.2) +
  stat_summary(fun.y=mean, geom="point")+
  scale_y_continuous(limits = c(0,1.2), breaks = seq(0,1.2,by= 0.3))+ theme_bw()+
  facet_wrap(.~protein)

ggsave("fractionation_westernblot.pdf", width = 5, height = 3, units="in")

sicon_values <- subset(inputD, status=="siCon" & protein=="HDAC2" & fraction == "N")$value
sitxnip_values <- subset(inputD, status=="siTXNIP" & protein=="HDAC2" & fraction == "N")$value
ttest <- t.test(sicon_values, sitxnip_values, alternative="two.sided", paired=T)
ttest

## paired two-sided T-test
# Cytoplasm : txnip : 0.02209 / hdac2 : 0.7904
# Nucleus : txnip : 0.03169 / hdac2 : 0.3471
## 
