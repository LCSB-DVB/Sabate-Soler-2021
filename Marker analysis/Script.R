library(data.table)
library(readxl)
library(tidyverse)
library(janitor)
library(ggplot2)

setwd("P:/MG/Paper/IF_Stainings_ImageAnalysis/20200429 Rstudio")
getwd()

####################
##LOADING OF DATA ##
####################

ABC1P1 =read.table('P:/MG/Paper/IF_Stainings_ImageAnalysis/20200429 Rstudio/ObjectsC1P1.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
ABC1P2 =read.table('P:/MG/Paper/IF_Stainings_ImageAnalysis/20200429 Rstudio/ObjectsC1P2.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
ABC1P3 =read.table('P:/MG/Paper/IF_Stainings_ImageAnalysis/20200429 Rstudio/ObjectsC1P3.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
ABC1P4 =read.table('P:/MG/Paper/IF_Stainings_ImageAnalysis/20200429 Rstudio/ObjectsC1P4.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
ABC1P5 =read.table('P:/MG/Paper/IF_Stainings_ImageAnalysis/20200429 Rstudio/ObjectsC1P5.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
ABC1P3eA =read.table('P:/MG/Paper/IF_Stainings_ImageAnalysis/20200429 Rstudio/ObjectsC1P3eA.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)
ABC1P3eB =read.table('P:/MG/Paper/IF_Stainings_ImageAnalysis/20200429 Rstudio/ObjectsC1P3eB.csv',header=TRUE, sep=',',stringsAsFactors = FALSE)

#########################
## MANIPULATING TABLES ##
#########################
#Adding columns
#1
ABC1P1$Batch = 1
ABC1P2$Batch = 2
ABC1P3$Batch = 3
ABC1P4$Batch = 4
ABC1P5$Batch = 5
ABC1P3eA$Batch= 3
ABC1P3eB$Batch = 3

ABC1list = list()
list.files(pattern = "*.csv") %>%
  map(~fread(.)) -> ABC1list
#fread = file read, it opens the files
#map function creates a list, it's applying the fread function to all the elements in the list

Table_ABC1 = bind_rows(ABC1list)

#
x =list()
for (i in seq(ABC1list)){
  x[[i]] = mutate(ABC1list[[i]],Batch =i)
}