##### ssGSEA metabolic gene sets #####
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
# 使用subset函数选择给定的列。
df <- subset(df,select= rownames(group))
df<-df[rowSums(df[,-1])!=0,]
# 导入代谢通路相关的基因集
Metab_gene_set<-import("Data/ssGSEA-Metabolic-Gene-Sets.xlsx")
geneset<-split(Metab_gene_set$Metagene,Metab_gene_set$GeneSet)

#### 免疫细胞丰度计算
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

normalization <- function(x){return((x-min(x))/(max(x)-min(x)))}  
nor_gsva_matrix <- normalization(scale_gsva_matrix)

library(pheatmap)
# 设置宽度和高度
pdf("Result/ssGSEA-Metabolic-genesets-heatmap.pdf",width=12,height=12)
# 绘制热图
pheatmap(nor_gsva_matrix,
         show_colnames =F,
         cluster_rows = F,cluster_cols = F,
         annotation_col = group,
         cellwidth = 10, cellheight = 10,
         fontsize = 10,
         #gaps_row = c(12,20)
         )
dev.off()

