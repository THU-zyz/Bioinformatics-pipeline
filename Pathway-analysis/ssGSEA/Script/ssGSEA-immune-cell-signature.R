##### ssGSEA immune cell signature analysis #####
library(tidyverse)
library(GSVA)
library(rio)

setwd("~/Desktop/ssGSEA/")
rm(list=ls())
source("Script/ssGSEA-function-repo.R")

# 导入表达矩阵
exp<-import("Data/count.all.txt")
df<-read4ENSEMBL2SYMBOL(exp)

# 导入分组信息矩阵，设置好样本对应状态
group<-import("Data/Sample_information.csv")
row.names(group)<-group$sample
group<-group%>%select(group)

# 导入免疫细胞基因集
# 免疫基因集来自于文章PMID:28052254
# 分析方法来自于文章PMID:30560866
# 包含28种免疫细胞的基因集
ImmuneGeneSet<-import("Data/ssGSEA-PMID28052254-Signature-28.xlsx")
MetabGeneSet<-import("Data/ssGSEA-gene-sets.xlsx")
geneset<-split(ImmuneGeneSet$Metagene,ImmuneGeneSet$`Cell type`)

#### 免疫细胞丰度计算
# gsva 接受的参数格式：
# 基因集输入格式是list;
# 基因表达量输入格式是矩阵，不是data.frame - 行为基因/列为样本;
# 默认情况，参数kcdf = "Gaussian"，适用于对数转换的microarray、RNA-seq的log-CPMs、log-RPKMs或log-TPMs；
# 当输入的表达矩阵是RNA-seq的raw count时，这个参数应该设置为kcdf="Poisson"；
# 参数method指定用于估计每个样本的基因富集分数的方法，默认情况下method="gsva",因此需要重新调整为method="ssgsea"。
ssGSEA_matrix<-gsva(expr=as.matrix(df),
                    gset.idx.list = geneset,
                    method = 'ssgsea',kcdf = 'Poisson', abs.ranking=TRUE)

#### 结果可视化
# 对数据进行z-score处理，scale原本是按列进行归一化，在这里经过转置后实现对行进行归一化。
# 即对每一种细胞类型在不同样本之间进行归一化处理。
scale_gsva_matrix <- t(scale(t(ssGSEA_matrix)))
# 修剪最大值最小值范围，将最小值限定在-2；最大值限定在2。
scale_gsva_matrix[scale_gsva_matrix < -2] <- -2
scale_gsva_matrix[scale_gsva_matrix > 2] <- 2

## 重新排序,按照抑癌、促癌、其它的顺序展示不同细胞的情况。
anti_tumor <- c('Activated CD4 T cell', 'Activated CD8 T cell', 'Central memory CD4 T cell', 
                'Central memory CD8 T cell', 'Effector memeory CD4 T cell', 'Effector memeory CD8 T cell',
                'Type 1 T helper cell', 'Type 17 T helper cell', 'Activated dendritic cell',
                'CD56bright natural killer cell', 'Natural killer cell', 'Natural killer T cell')

pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'CD56dim natural killer cell',
               'Immature dendritic cell', 'Macrophage', 'MDSC', 'Neutrophil', 'Plasmacytoid dendritic cell')

anti <- rownames(scale_gsva_matrix) %in% anti_tumor 
pro <- rownames(scale_gsva_matrix) %in% pro_tumor
non <- !(anti|pro)
# 重新排列一下展示顺序,依次展现抗癌相关细胞、促癌相关细胞、以及其他细胞。
scale_gsva_matrix <- rbind(scale_gsva_matrix[anti,],scale_gsva_matrix[pro,],scale_gsva_matrix[non,])
normalization <- function(x){return((x-min(x))/(max(x)-min(x)))}  
nor_gsva_matrix <- normalization(scale_gsva_matrix)
# 使用subset函数选择给定的列。
nor_gsva_matrix <- subset(nor_gsva_matrix,select= rownames(group))

library(pheatmap)

pheatmap(nor_gsva_matrix,
         show_colnames =F,
         cluster_rows = F,cluster_cols = F,
         annotation_col = group,
         cellwidth = 10, cellheight = 10,
         fontsize = 10,
         gaps_row = c(12,20))




