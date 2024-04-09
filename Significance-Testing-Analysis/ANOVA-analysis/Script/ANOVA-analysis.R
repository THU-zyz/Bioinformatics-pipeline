##### ANOVA analysis #####
library(tidyverse)
library(xlsx)
# 如果是第一次打开只运行一次，报错则关闭Rstudio再重新打开此文件重新运行
rm(list = ls())

source("../Script/ANOVA-analysis-function-repo.R")
# 输入文件所在地址
file_path="Data/P7-PLSDA.xlsx"
# 输出结果文件名
output_name="P7"

# 输入数据需要满足以下条件,从第三列开始是记录了代谢物丰度的列，前两列分别是样本名称和分组group
ANOVA_analysis(file_path,output_name)
