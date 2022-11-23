## Written by KyungTae Lee
## Last update : 2017-09-17
## R script : "ROC_average_curves.R"
## Script to plot average ROC cuvers using positive and negative PPI sets.

setwd("Z:/dataset/UW_protein/dataSummary_result/2017_03_08_Drosophila_Human/Human_protein/SAINT_analysis/evaluation/result/top3Controls")
getwd()

prefix<- 'Human_' ## output prefix
cutoff_type <- "saintScore" # Type of cutoffs to be used to plot the ROC cuvers
if(!require(ROCR)){
    install.packages("ROCR")
}
library(ROCR)
library(RColorBrewer)

true_table<- read.table("true_predictions.SAINTscore.txt", header= TRUE, sep= "\t")
## true_table should have "Prediction" as header and prediction value at the first column
true_table <- data.frame(true_table)

false_table <-read.table("false_predictions.1000RandomSet.SAINTscore.txt",header= TRUE, 
                         sep= "\t")
## false_talbe should have multiple columns of prediction values that are randomly sampled. Length of each column should equal to that of true_table

false_table <- data.frame(false_table)
head(true_table)
head(false_table)

if (FALSE){
  get_negative_preds <- function(set_size, input_negative_preds){
    random_negative_pred <- sample(input_negative_preds, set_size)
    return (random_negative_pred)
    }
}

true_preds <- true_table$prediction
true_num<- length(true_preds)
true_num

true_labels<- rep(c(1), true_num)
false_labels <- rep(c(0), true_num)
class_labels<- c(true_labels, false_labels)
## Making class label vector. "1"= TRUE, "0"= FALSE

multiple_class_labels<- list()
multiple_preds <-list()
## Initialize lists that will be used as inputs for ROCR package

for (i in 1:1000){
  multiple_class_labels[[i]]<- class_labels
  negative_pred <- false_table[[i]]
  pred_set <- c(true_preds, negative_pred)
  multiple_preds[[i]]<- pred_set  
}

## Making inputs for ROCR (multiple sets, cross-validation set for 
## example)
pred<- prediction(multiple_preds, multiple_class_labels)
roc.perf= performance(pred,"tpr", "fpr")
attributes(roc.perf)$alpha.values 

##Code to set maximum threshold to 1
if (FALSE){
  if (cutoff_type=="SpecThreshold" | 
      cutoff_type=="SpectralCount") {
  } else {
    for (i in 1:1000){
      attributes(roc.perf)$alpha.values[[i]][1]<-1
    }  
  }
}

## To avoid error of cutoffs biggern than 1.0 appearing. Crucial step!

if (FALSE){
  plot(roc.perf, col="grey82", lwd = 2,
       main = paste("ROC curves (1000 sampling).",prefix, sep= ""),
       cex.lab=1.5)
  
  #cex.lab=1.5, cex.axis= 1.5)
  axis(2, cex.axis=1.5)
  axis(1, cex.axis=1.5)
}


pdf(paste(prefix,".multiple_roc.pdf", sep=""))
plot(roc.perf, col="grey82", lwd = 2,
     main = paste("ROC curves (1000 sampling).",prefix, sep= "" ),
     cex.lab= 1.5)
plot(roc.perf,lwd=3,avg="threshold",spread.estimate="stddev",add=TRUE)
axis(2, cex.axis=1.5)
axis(1, cex.axis=1.5)
abline(a=0, b=1)
dev.off()

summary(roc.perf)

averageRoc <- function(prefix, perf_object) {
  color_pal= rev(rainbow(256, start = 0, end = 4/6))
  
  pdf(paste(prefix,"average_roc","pdf",
            sep="."))
  plot(perf_object, lwd=3, avg="threshold", 
       spread.estimate= "stddev",
       colorize.palette= color_pal,
       colorize= TRUE,  colorkey.relwidth=0.25,
       colorkey.pos= "right",
#      print.cutoffs.at=c(3.75,-1.12586),   
       # 5'ss optional cutoff
       print.cutoffs.at=c(0.85),
       # 3'ss optional cutoff

       main= paste("Average ROC (1000_sampling).", 
                   prefix, sep="."),
       cex.lab=1.5)
  axis(2, cex.axis=1.5)
  axis(1, cex.axis=1.5)
  roc.auc<- performance(pred, measure= "auc")

  auc_values<-unlist(roc.auc@y.values, use.names = FALSE)
 
  auc_mean<- mean(auc_values)
  auc_mean <- round(auc_mean, digits = 3)
  auc_mean_legend<- paste("Mean AUC = ", 
                          toString(auc_mean), sep="")
  text(0.7,0.3,labels= auc_mean_legend, cex=1)
  abline(a=0, b=1)
  dev.off()
}

averageRoc(prefix, roc.perf)

opt.cut= function(perf,pred){
  cut.ind= mapply(FUN=function(x,y,p){
    d=(x-0)^2 + (y-1)^2
    ind= which(d== min(d))[1]
    c(sensitivity= y[[ind]], specificity= 1-x[[ind]],
      cutoff= p[[ind]])
  }, perf@x.values, perf@y.values, pred@cutoffs)
}  

print (opt.cut(roc.perf, pred))

optimal_cutoffs<- opt.cut(roc.perf, pred)
cutoff<- mean(optimal_cutoffs[3,])
print (cutoff)
#print (paste("Average Cutoff", toString(cutoff), sep = " : "))
#write.table(optimal_cutoffs, "test.txt", sep = "\t", row.names = T)





