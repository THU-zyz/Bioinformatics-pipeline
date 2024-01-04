library(tidyverse)
library(parallel)
library(GenomicFeatures)
# 设置启用计算机3/4的核
cl<-makeCluster(0.75*detectCores())
# BiocManager::install("GenomicFeatures")
# setwd("~/Desktop/博士生活/科研生活/Hu_lab/Github/Bioinformatics-pipline/scripts/count2TPM/")

# 1) 从featurecounts中得到有效基因长度
expr_df<-read.table("../Data/A-2333.featurecounts.all.txt",header=T,check.names=F)
gene_efflen<-expr_df%>%dplyr::select(ENSEMBL=Geneid,efflen=Length)
write.csv(gene_efflen,"gene_len.csv",row.names = F)

# 2)从gtf文件中计算有效基因长度
hg38db<-makeTxDbFromGFF("../Data/Homo_sapiens.GRCh38.110.gtf",format="gtf")
# 提取外显子
exons_gene<-exonsBy(hg38db,by="gene")
head(exons_gene)
# 计算总外显子长度：用reduce去掉重复冗余的部分，width统计长度，最后计算总长度。
exons_gene_lens<-parLapply(cl,exons_gene,function(x){sum(width(reduce(x)))})
exons_gene_lens[1:10]
# 转换为dataframe
gene_efflen_Homo_sapiens.GRCh38.110<-data.frame(EMSEMBL=names(exons_gene_lens),
                          efflen=as.numeric(exons_gene_lens))
# 对比了一下这两种方法拿到的efflen结果一样
write.csv(gene_efflen_Homo_sapiens.GRCh38.110,"gene_efflen_Homo_sapiens.GRCh38.110.csv",row.names=F)

####### convert raw reads counts to cpm/tpm/fpkm
# 注意gene_matrix需要满足一下格式，第一列为ENSEMBL,记录基因的ENSEMBL id；
# 后续列为每个样本在这个基因上的read counts数
# define count2TPM function
count2tpm<-function(gene_matrix,gene_efflen){
  gene_matrix<-gene_matrix%>%inner_join(gene_efflen,by="ENSEMBL")
  # 将非表达列
  gene_matrix<-gene_matrix%>%dplyr::select(ENSEMBL,efflen,everything())
  expr_df<-gene_matrix[,-(1:2)]
  # convert the data type of the expr_df columns to numeric.
  expr_df<-apply(expr_df,2,as.numeric)
  # 基因长度，目标基因的外显子长度之和除以1000，单位是kb，不是bp
  kb<-gene_matrix$efflen/1000
  # per thousand scaling factor 每千碱基reads 长度标准化。
  rpk<-expr_df/kb
  # per million scaling factor 每百万缩放因子 
  tpm <-t(t(rpk)/colSums(rpk)*1000000)
  gene_matrix[,-(1:2)]<-tpm
  gene_matrix<-gene_matrix%>%dplyr::select(-efflen)
  return(gene_matrix)
}

# define count2fpkm function
count2fpkm<-function(gene_matrix,gene_efflen){
  gene_matrix<-gene_matrix%>%inner_join(gene_efflen,by="ENSEMBL")
  gene_matrix<-gene_matrix%>%dplyr::select(ENSEMBL,,efflen,everything())
  expr_df<-gene_matrix[,-(1:2)]
  # convert the data type of the expr_df columns to numeric.
  expr_df<-apply(expr_df,2,as.numeric)
  # 基因长度，目标基因的外显子长度之和除以1000，单位是kb，不是bp
  kb<-gene_matrix$efflen/1000
  # per thousand scaling factor 每千碱基reads 长度标准化。
  rpk<-expr_df/kb
  fpkm<-t(t(rpk)/colSums(expr_df)*10^6)
  gene_matrix[,-(1:2)]<-fpkm
  return(gene_matrix)
}

# define count2cpm function
count2cpm<-function(gene_matrix){
  expr_df<-gene_matrix[,-1]
  # convert the data type of the expr_df columns to numeric.
  expr_df<-apply(expr_df,2,as.numeric)
  cpm<-t(t(expr_df)/colSums(expr_df)*1000000)
  gene_matrix[,-1]<-cpm
  return(gene_matrix)
}
