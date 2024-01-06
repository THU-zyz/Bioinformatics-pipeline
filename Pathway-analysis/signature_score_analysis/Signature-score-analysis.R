library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(ggsci)

rm(list=ls())
setwd("./")
# Invoke custom fuction library
source('Script/count2TPM.R')
source('Script/signature_score_function_repo.R')

##### parameter settings #####
# File storage path to expression matrix
expression_path="Data/count.all.txt"
# File storage path to pathway gene
pathway_path="Data/Embryonic-Pausing-Signature.xlsx"
# File storage path to sample information
sample_info_path="Data/Sample_information.csv"
# File storage path to gene effective length
gene_efflen_path<-"Data/gene_efflen.csv"
# custom color palette
my_palette=c("T1"="#000000","T2"="#0000CD","T3"="#B22222","T5"="#DAA520","T6"="#228B22")
# comparisons<-list(c("T1","T2"),c("T1","T3"),c("T1","T5"),c("T1","T6"))
# Significance testing of group information.
comparisons=".all."
# pathway name 
pathway_name="Embryonic-Pausing-Signature"
pathway_title="Embryonic Pausing Signature"
# Dataset name
data_name="Embryo"


##### Read expression matrix and signature gene list.
##### Then, carry out data preprocessing
run_signature_score_analysis(expression_path,pathway_path,sample_info_path,gene_efflen_path,my_palette,comparisons,pathway_name,pathway_title,data_name)





