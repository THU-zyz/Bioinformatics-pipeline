##### sgRNA integrative analysis pipeline
# 安装R包
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("MAGeCKFlute")

# 导入包
library(MAGeCKFlute)
library(tidyverse)
library(ggplot2)
library(ggrepel)
rm(list=ls())
# 绝对路径
# setwd("~/Desktop/sgRNA-analysis/Mageck-Test-analysis/")
# 相对路径, 打开这个脚本时所在路径的相对路径。
setwd("../../")
source("Script/sgRNA-Post-Processing/sgRNA-Integrative-Analysis-Function-Warehouse.R")
# 设置sgrna或gene summary所在路径。
gene_summary_path = "Data/sqh20231106/WT-NK/Mageck-test/WT_vs_NK.gene_summary.txt"
sgrna_summary_path = "Data/sqh20231106/WT-NK/Mageck-test/WT_vs_NK.sgrna_summary.txt"

##### gene summary 结果文件数据分析 #####
# 基因列表，用于在图中展示这些基因。
gene_list = c("Pvr","Nceh1","Tfrc","Slc31a1","Ppic","Neo1","Tmed2","Stt3a")
# 标题名称。
title_name = 'WT vs NK' 
output_path = "Result/sqh20231106/sgRNA-Integrative-Analysis/"
# 当LFC_threshold>0时，选择正向筛选的结果；
# 当LFC_threshold<0时，展示负向筛选的结果；
LFC_threshold = -2
# 设置自定义颜色盘,与基因list对应。
my_color_paltte = c("Nceh1" = "#FA7F6F","Pvr"="#C76DA2","Tfrc"="#FFBE7A","Slc31a1"="#BEB8DC","Ppic"="#FFFF77","Neo1"="#8ECFC9","Tmed2"="#99FFFF","Stt3a"="#82B0D2")
#  检查并创造结果文件夹
create_output_dir4sgRNA(output_path)

sgRNA_analysis4gene(gene_summary_path, output_path, gene_list, my_color_paltte, LFC_threshold, title_name)

##### sgrna summary 结果文件数据分析 ##### 
gene_summary_path = "Data/sqh20231106/WT-NK/Mageck-test/WT_vs_NK.gene_summary.txt"
sgrna_summary_path = "Data/sqh20231106/WT-NK/Mageck-test/WT_vs_NK.sgrna_summary.txt"
output_path = "Result/sqh20231106/sgRNA-Integrative-Analysis/"
title_name = 'WT vs NK' 
binwidth = 0.1
interval = 0
sgRNA_analysis4sgrna(sgrna_summary_path,gene_summary_path, gene_list, output_path, title_name, binwidth = 0.1,interval = 0)
