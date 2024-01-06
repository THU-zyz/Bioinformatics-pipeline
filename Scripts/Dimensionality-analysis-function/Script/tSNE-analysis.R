#####  t-SNE analysis #####
# t-dsitributed Stochastic Neighbor Embedding. 
library(tidyverse)
library(Rtsne)
library(openxlsx)
rm(list=ls())
setwd("../")
source("Script/Dimensionality-analysis-function-repo.R")
##### Parameters settings #####

files<-"Data/PDX-KW-test.xlsx"
sheetnames="Sheet1"
output_name<-"PDX-t-SNE-Result"
my_palatte<-c("T1" = "#6fe7dd","T2" = "#3490de", "T3" = "#6639a6","T4"="#FF8C00","T5"="#6B8E23","T6"="#A0522D")

# run tSNE_analysis
tSNE_analysis(files,sheetnames,output_name,my_palatte,perplexity=5,normalized=T)

