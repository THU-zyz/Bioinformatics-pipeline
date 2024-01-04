library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)
library(ggsci)

rm(list=ls())
setwd("~/Desktop/DTP/signature_score_analysis/")
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
pathway_name="EmbryonicPausingSignature"
plot_title="Embryonic Pausing Signature"
# Dataset name
data_name="Embryo"


##### Read expression matrix and signature gene list.
##### Then, carry out data preprocessing
## 1) Data loading section.
# 1.1) gene expression matrix
expression<-read.table(expression_path,header = TRUE)
expression<-expression%>%dplyr::select(ENSEMBL=gene_id,everything())

# 1.2) To establish a general "embryonic pausing signature",
# First measured DE in each embryonic model independently
# Gene sigificantly dysregulated in both models were selected for the signature
# Other Gene signatures such as mTORC1 response were obtained from the MSigDB "Hallmark gene sets"
sign_gene<-read.xlsx(pathway_path,sheet = 1)
sign_gene<-sign_gene%>%dplyr::select(SYMBOL=Gene,everything())
# the expression value was multiplied by -1 if the gene was downregulated in the embryonic models,
# while the value of genes upregulated in the embryonic models was kept unchanged，
### ==so that a positive expression value would always reflect a change similar to embryonic pausing.== 

# 1.3) create a sample imformation matrix
sample_info<-read.csv(sample_info_path,header=T,)

# 1.4) Read the matrix of effective gene lengths for calculating TPM
gene_efflen<-read.csv(gene_efflen_path)

## 2) Data processing section.
# 2.1) perform data processing on gene expression matrix.
gene_expression<-gene_matrix_preprocessing(expression,sign_gene,sample_info,gene_efflen)
# 2.2) calculation signature score.
sign_score<-signature_score(gene_expression,sample_info)
# save sign_score
write.csv(sign_score,paste0("Result/Signature_score_for_",pathway_name,"_using_",data_name,".csv"),row.names = F)
# 2.3) plot signature boxplot
plot1<-signature_score_plot(sign_score,title=plot_title,my_palette,comparisons=comparisons,ylab="Signature Score",Test_Methods="t.test",hide_ns=T,sign_label="p.signif")
# 2.4) save the figure.

ggsave(plot1,device = cairo_pdf,
       path="Result/Figure/",
       filename = paste0(data_name,"_signature_score_for",pathway_name,"_plot.pdf"),width = 8,height = 10,units ="in")





