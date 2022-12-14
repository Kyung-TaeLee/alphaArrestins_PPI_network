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

## 2. Ngoc fractionation western blot

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

## 3. Ngoc RT-qPCR and RNA-seq average values

#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

inputD <- read.delim("qpcr_tubulin.RNAseq.replicates.txt", header=T, sep="\t", check.names = F)
head(inputD)

avger_inputD <- data_summary(inputD, "value", groupnames= c("condition","method"))

# Use position=position_dodge()
p<- ggplot(data=avger_inputD, aes(x=condition, y=value)) +
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd), width=.2,
                position=position_dodge(.9)) + 
  theme_bw()+
  facet_wrap(.~method, scales = "free_y")+
  geom_point(data= inputD, aes(x=condition, y= value))
p

## RT-qPCR t test (two-sided, paired)
sicon_values <- subset(inputD, condition=="siCon" & method=="rt-qpcr")$value
sitxnip_values <- subset(inputD, condition=="siTXNIP" & method=="rt-qpcr")$value
ttest <- t.test(sicon_values, sitxnip_values, alternative="two.sided", paired=T)
ttest
# p-value = 0.0001062

ggsave("qpcr_tubulin.RNAseq.replicates.pdf")
