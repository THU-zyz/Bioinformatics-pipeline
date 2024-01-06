#####  PLS-DA analysis #####
# Principal Component Analysis PCA analysis.
library(tidyverse)
library(Rtsne)
library(openxlsx)
library(mixOmics)
rm(list=ls())
setwd("../")
source("Script/Dimensionality-analysis-function-repo.R")

##### Parameters settings #####
files<-"Data/PDX-KW-test.xlsx"
sheetnames="Sheet1"
output_name<-"PDX-PLS-DA-Result"

# run plsda_analysis
plsda_analysis(files,sheetnames,output_name)



