##### sgRNA library quality assessment #####
library(tidyverse)
library(openxlsx)
rm(list=ls())
# 绝对路径
# setwd("~/Desktop/sgRNA-analysis/Mageck-Test-analysis/")
# 相对路径
setwd("../../")
source("Script/sgRNA-Post-Processing/sgRNA-Data-Processing-Function-Warehouse.R")
##### parameters  settings #####
input_files="Data/jxj20240110/mageck_count/PCR2.count.txt"
output_dir="Result/jxj20240110/"
plot_name="Log2 Readcount"
sample_info_path="Data/jxj20240110/PCR2-sample-info.xlsx"

# sgRNA结果绘制Log2Readcount矩形图。
sgRNA_Log2Readcount(input_files,output_dir,sample_info_path,plot_name)

# 读取的count.txt文件中包含多个组的样本
input_files="Data/wms20231023/Plasmid/Mageck-count/plasmid.count.txt"
output_dir="Result/wms20231023/Plasmid/"
plot_name="Log2 Readcount"
sample_info_path="Data/wms20231023/Plasmid/Plasmid-sample-info.xlsx"

sgRNA_Log2Readcount(input_files,output_dir,sample_info_path,plot_name)

# 读取的count.txt文件中包含多个组的样本
input_files="Data/wms20231023/PC9/Mageck-count/PC9_all.count.txt"
output_dir="Result/wms20231023/PC9/"
plot_name="Log2 Readcount"
sample_info_path="Data/wms20231023/PC9/PC9-sample-info.xlsx"

sgRNA_Log2Readcount(input_files,output_dir,sample_info_path,plot_name)

