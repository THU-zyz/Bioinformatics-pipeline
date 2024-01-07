##### Fuzzy Clustering of Time Series Gene Expression Data #####
##### Fuzzy C-Means Clustering FCM #####
# Analyze the temporal trends of gene or protein expression in transcriptomic, proteomic, and metabolomic data with time series characteristics.
# Group genes, proteins, and metabolites with similar expression patterns into one category to help understand the dynamic patterns of these biological molecules and their association with function.
# Mfuzz requires the data provided to be normalized expression levels, such as TPM or FPKM. When there are biological replicates, and a single time point corresponds to multiple samples, it is necessary to calculate the average expression level for each group.

library(tidyverse)
library(Mfuzz)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
rm(list=ls())
setwd("Mfuzz-analysis/")
source("Script/Mfuzz-analysis-function-repo.R")

# For metabolomic data, it can be directly imported for calculating
# For the raw data of transcriptomics, the transcriptomic data needs to be converted into TPM format before calculation. 

# Before Mfuzz analysis, a preliminary screening is needed to select a subset from all genes, metabolites, and proteins.
# The, perform cluster analysis on this subset. The selection criterion is P<0.05 after KW test (for metabolomic data, choose P<0.1).

##### Parameters Settings #####
# Type.I Metabolomics data
# For metabolomic data
paras<-list()
paras[["type"]]<-"Metabolomics"
files="Data/Metabolomics/PDX-KW-test.xlsx"
sheetnames<-"Sheet1"
output_name<-"PDX-Mfuzz-Meta"
pvalue_threshold<-0.1
Mfuzz_number<-4
labels<-c("Untreated","Phase1","Phase2","Phase3","Phase4","Phase5")


Metabolomics_Mfuzz<-Mfuzz_KW(paras,files,sheetnames,output_name,pvalue_threshold,labels,Mfuzz_number)


# Type.II Transcriptions data
# For the raw data of transcriptomics
paras<-list()
paras[["type"]]<-"Transcriptomics"
paras[["sample_info_path"]]<-"Data/Transcriptomics/Sample_information.csv"
paras[["gene_efflen_path"]]<-"Data/Transcriptomics/gene_efflen.csv"
files="Data/Transcriptomics/PDX-reads-count.txt"
sheetnames="Sheet1"
output_name<-"PDX-Mfuzz-RNA"
pvalue_threshold<-0.05
Mfuzz_number<-4
labels<-c("Untreated","Phase1","Phase2","Phase4","Phase5")

Transcriptions_Mfuzz<-Mfuzz_KW(paras,files,sheetnames,output_name,pvalue_threshold,labels,Mfuzz_number)


