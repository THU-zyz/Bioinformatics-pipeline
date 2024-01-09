##### RNA-seq differential expression analysis #####
library(DESeq2)
library(tidyverse)
library(ggplot2)
library(openxlsx)
# Human Genome Annotation
library(org.Hs.eg.db)
# Mouse genome Annotation
library(org.Mm.eg.db)
library(clusterProfiler)
library(ggrepel)
rm(list=ls())
# setting the script's execution directory
setwd("../")
source("Script/RNA-seq-DEA-function-repos.R")


##### Part.I The data contains only one group or requires differential gene analysis for just a specified single group #####
# 1.1 Parameter settings 
# Data input address.
input_files="Data/count.all.txt"
# Sample information address. 
sample_info_path="Data/Sample_information.csv"
# Set the two groups to be compared.
comparisons_T12<-c("T1","T2")
# Set one group as the control group.
refgroup="T1"
# Set the output name
output_name="T12"
# Human metabolic genes address.
human_metabolic_gene="Data/human_metabolic_genes.csv"
# Set the title of the volcano plot.
plot_title="Phase1 vs Untreated"
# Set the title of the volcano plot for xlab.
xlab="log2(FoldChange,Phase1/Untreated)"
# Set the title of the volcano plot for ylab.
ylab="-log10(FDR)"
# Set the FDR threshold value.
padjThread=0.05
# Set the threshold for the total number of reads
readssumThread=1
# Set the log2Foldchange threshold value.
L2FCThread=1.5


# 1.2 Perform differential expression analysis.
allNmr<-Differential_expression_analysis(input_files,sample_info_path,
      comparisons_T12,refgroup,output_name,human_metabolic_gene,plot_title,xlab,ylab,
      padjThread=0.05,readssumThread=1,L2FCThread=1.5)
# Separately extract all genes and those related to metabolism.
all_gene<-allNmr$all
mr_gene<-allNmr$mr
# filter out differential expression gene.
all_DE_gene<-all_gene%>%filter(sig!="none")
mr_DE_gene<-mr_gene%>%filter(sig!="none")

# Perform KEGG enrichment for all genes.
all_kegg_enrich<-KEGGpathwayEnirch(gene_set=all_DE_gene$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db",organism="hsa",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
# Perform KEGG enrichment for metabolism-related genes
mr_kegg_enrich<-KEGGpathwayEnirch(gene_set=mr_DE_gene$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db",organism="hsa",pvalueCutoff = 0.05,qvalueCutoff = 0.05)

# Draw a bubble chart containing all the pathways enriched for all gene.
KEGG_enrich_plot(all_kegg_enrich,paste0(output_name,"-all-"))
# Draw a bubble chart containing metabolism-related the pathways enriched for all gene.
KEGG_enrich_metabolism_related_plot(all_kegg_enrich,paste0(output_name,"-all-"),pathway_list_path)
# Draw a bubble chart containing all the pathways enriched for metabolism-related gene.
KEGG_enrich_plot(mr_kegg_enrich,paste0(output_name,"-mr-"))
# Draw a bubble chart containing metabolism-related the pathways enriched for metabolism-related gene.
KEGG_enrich_metabolism_related_plot(mr_kegg_enrich,paste0(output_name,"-mr-"),pathway_list_path)


#####  Part.II Draw a KEGG enrichment chart containing multiple groups. #####
# To compare the KEGG enrichment results of multiple groups, they need to be integrated and then drawn in a bubble chart format with the group names as row names.
# Data input: Requires a dataframe that includes a column named 'cluster', recording the grouping information of the integrated KEGG enrichment results.
# 2.1 Parameter settings.
# Data input address.
input_files="Data/count.all.txt"
# Set the sample information file.
sample_info_path="Data/Sample_information.csv"
# Set the two groups to be compared.
comparisons_dir=list(T12=c("T1","T2"),T13=c("T1","T3"),T15=c("T1","T5"),T16=c("T1","T6"))
# Set which group to be designated as the control group.
refgroup="T1"
# Set the output file name.
outputname_dir=list(T12="T12",T13="T13",T15="T15",T16="T16")
# Set the name of the graph.
plot_name="DTP"
# Set the address of the human metabolic gene file.
human_metabolic_gene="Data/human_metabolic_genes.csv"
# Set the title of the volcano plot.
plot_title=""
# Set the title of the volcano plot's x-axis.
xlab="log2(FoldChange)"
# Set the title of the volcano plot's y-axis.
ylab="-log10(FDR)"
# Run the multi group KEGG pathway enrichment function.
multi_group_KEGGpathwayEnrich(input_files,sample_info_path,
                                        comparisons_dir,refgroup,outputname_dir,plot_name,
                                        human_metabolic_gene,plot_title,xlab,ylab,
                                        padjThread=0.05,readssumThread=1,L2FCThread=1.5)



