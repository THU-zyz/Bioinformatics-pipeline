####### Kruskal-Wallis test，a non-parametric test for comparing variances among multiple independent samples. ######
# 检查并安装缺失的R包
packages_needed=c("readxl","tidyverse","dunn.test")
for (package in packages_needed) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package)
  }
}

# import R packages for analysis 
library(readxl)
library(tidyverse)
library(dunn.test)
setwd("../") 
# Strings are not automatically converted into factors by default.
options(stringsAsFactors = F)
source('Script/KW—test-function-repo.R')
##### parameters settings #####
# Setting the file name for reading.
files="Data/Metabolomics.xlsx"
# Setting the sheet name for reading.
sheetnames="Sheet1"
# Setting the output file name.
output="Result/embryo_KW_test_result.csv"

### Run KW test function.
Result<-KW_test(files = files,sheetnames = sheetnames,output = output)
