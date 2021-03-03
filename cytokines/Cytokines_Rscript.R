---
#title: "Metabolites_secretion"
#author: "Sònia"
#date: "2020 M06 26"
#output: html_document
---

{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)




setwd("P:/MG/Paper/CK pannel")
metadata <- read.csv("Cytokines.csv",header=T, row.names=1, sep=",")




#install.packages("matrixStats")
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(matrixStats)
```


##No normalization of dataset
#Just to see whether replicates group together
pheatmap(as.matrix(metadata[,1:12]), color = colorRampPalette(rev(brewer.pal(n = 12, name = "Spectral")))(100),legend=TRUE, main="", angle_col = 45, fontsize_row = 8, fontsize_col=10,cellwidth = 15, width=15, height = 10)


#Transform Row index into 1st column called "Metabolite"
  library(dplyr)
data <- tibble::rownames_to_column(metadata, "Cytokine")




#With the pooled data for a pooled heatmap
data1<-data%>%
  mutate(Midbrain_organoids = rowMedians(as.matrix(data[,2:4])))%>%
  mutate(Assembloids_K7 = rowMedians(as.matrix(data[,5:7])))%>%
  mutate(Assembloids_163 = rowMedians(as.matrix(data[,8:10])))%>%
  mutate(Assembloids_EPI = rowMedians(as.matrix(data[,11:13])))
  
data2<- column_to_rownames(data1,var="Cytokine")
data3<-data2[,13:16]


#rename cell lines
#colnames(data3)<-c('Ctrl_1', 'Ctrl_2', 'Ctrl_3', 'IPD_1', 'IPD_2', 'IPD_4', 'IPD_3', 'IPD_5', 'IPD_6')

#normalization
Z_Normalize<-function(x){return((x-mean(x))/(sd(x)))}

NormalisedData_Zscore = t(apply(as.matrix(data2),1,Z_Normalize))

pheatmap(as.matrix(NormalisedData_Zscore), color = colorRampPalette(rev(brewer.pal(n = 8, name = "Spectral")))(100),legend=TRUE, main="Cytokine release", angle_col = 45, fontsize_row = 7, fontsize_col=10,cellwidth = 30, width=30, height = 15)


#Now with pooled samples
Z_Normalize<-function(x){return((x-mean(x))/(sd(x)))}

NormalisedData_Zscore = t(apply(as.matrix(data3),1,Z_Normalize))

pheatmap(as.matrix(NormalisedData_Zscore), color = colorRampPalette(rev(brewer.pal(n = 8, name = "Spectral")))(100),legend=TRUE, angle_col = 45, fontsize_row = 11, fontsize_col=13,cellwidth = 30, width=30, height = 15)


#MOs vs pMGL-MOs
datapooled<-data%>%
  mutate(MO_wo = rowMedians(as.matrix(data[,2:4])))%>%
  mutate(pMGLMO = rowMedians(as.matrix(data[,5:13])))

datapooled1<- column_to_rownames(datapooled,var="Cytokine")
datapooled2<-datapooled1[,13:14]

Z_Normalize<-function(x){return((x-mean(x))/(sd(x)))}

NormalisedData_Zscore = t(apply(as.matrix(datapooled2),1,Z_Normalize))

pheatmap(as.matrix(NormalisedData_Zscore), color = colorRampPalette(rev(brewer.pal(n = 8, name = "Spectral")))(100),legend=TRUE, main="Cytokine release", angle_col = 45, fontsize_row = 7, fontsize_col=10,cellwidth = 30, width=30, height = 15)






