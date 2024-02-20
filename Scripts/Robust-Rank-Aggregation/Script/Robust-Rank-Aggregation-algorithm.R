##### Robust Rank Aggregation RRA
# RRA method is based on the assumption that each gene in each dataset is randomly arranged.
# If the gene ranks high in all datasets, the associated P-value is lower, the possibility of differential gene expression is greater.
# Through rank analysis, 343 integrated DEGs, consisting of  111 upregulated genes and 232 downregulated genes, were identified by the RRA method.
BiocManager::install("RobustRankAggreg",ask=F,update=F)
library(RobustRankAggreg)
library(clusterProfiler)
library(pheatmap)
library(openxlsx)
library(pheatmap)
library(org.Hs.eg.db)
# Method implementing various gene list aggregation methods, most notably Robust Rank Aggregation
# aggregateRanks函数对于多个排好序的基因集，在求交集的同时还考虑它们的排序情况；
# 挑选在多个数据集上都表现差异的基因，并每次差异都排名靠前的那些。
# For example
set.seed(10086)
ENSEMBL2SYMBOL<-function(gene_matrix,OrgDb = "org.Hs.eg.db"){
  gene_SYMBOL<-bitr(gene_matrix$ENSEMBL, fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = "org.Hs.eg.db")
  df<-gene_matrix%>%inner_join(gene_SYMBOL,by="ENSEMBL")
  df<-df%>%dplyr::select(SYMBOL,everything(),-ENSEMBL)
  return(df)
}

glist<-list(sample(letters,16),
            sample(letters,15),
            sample(letters,13))
freq = as.data.frame(table(unlist(glist)))
# Aggregate the inputs
ag = aggregateRanks(glist = glist,
                    N = length(letters))
ag$Freq = freq[match(ag$Name, freq$Var1),2]

# For reality gene list example
set.seed(12345678)
data(geneList,package="DOSE")
head(geneList)
deg=data.frame(gene=names(geneList),
               logFC=as.numeric(geneList),P.Vaule=0)

expression_status<-function(df){
  # 判断是否pvalue>0.05;如果满足则为stable基因
  df$padj<-as.numeric(df$padj)
  df$sig=ifelse(df$padj>0.05,'none',
              # 判断padj < 0.05的基因中是否满足log2FoldChange > 1.5,满足则为up(上调).
              ifelse(df$log2FoldChange>1.5,'up',
                     # 判断log2FoldChange < -1.5的基因,满足则为down(下调).
                     ifelse(df$log2FoldChange < -1.5,'down','none')))
  return(df)
}
get_sig<-function(df,sig_type){
  df<-df%>%filter(sig==sig_type)
  df<-df[order(df$log2FoldChange,decreasing = T),]
  return(df[df$sig==sig_type,'SYMBOL'])
}

glist=list(expression_status(df1),
           expression_status(df2),
           expression_status(df3))
glist_RRA<-aggregateRanks(glist=glist,N=length(unique(unlist(glist))))
tmp = as.data.frame(table(unlist(glist)))
glist_RRA$Freq=tmp[match(RRA$Name,tmp[,1]),2]

# plot heatmap
qualified_RRA_gene<-RRA[RRA$Score < 0.05,1]
qualified_RRA_genelist<-dataframe(
  deg2=deg2[qualified_RRA_gene,'logFC'],
  deg3=deg3[qualified_RRA_gene,'logFC'],
  deg4=deg4[qualified_RRA_gene,'logFC']
)
rownames(updat)=gs
pheatmap(updat,display_numbers = T)

setwd("~/Desktop/博士生活/科研生活/Hu_lab/Github/Bioinformatics-pipeline/Scripts/Robust-Rank-Aggregation/")
DEG1<-read.xlsx("Data/DEG1.xlsx")
DEG2<-read.xlsx("Data/DEG2.xlsx")
DEG3<-read.xlsx("Data/DEG3.xlsx")
glist<-list(DEG1=DEG1,
            DEG2=DEG2,
            DEG3=DEG3)

RRA_analysis<-function(glist){
  up_gene_list<-list()
  down_gene_list<-list()
  gene_list_log2FC<-list()
  for( key in names(glist)){
    item<-glist[[key]]
    item<-ENSEMBL2SYMBOL(item,OrgDb = "org.Hs.eg.db")
    item<-expression_status(item)
    up_gene<-get_sig(item,"up")
    down_gene<-get_sig(item,"down")
    up_gene_list[[key]]<-up_gene
    down_gene_list[[key]]<-down_gene
    gene_list_log2FC[[key]]<-item%>%dplyr::select(SYMBOL,log2FoldChange)
    
  }
  # 计算上调基因的RRA排名
  up_gene_RRA<-aggregateRanks(glist=up_gene_list,N=length(unique(unlist(up_gene_list))))
  up_gene_freq<-as.data.frame(table(unlist(up_gene_list)))
  up_gene_RRA$Freq<-up_gene_freq[match(up_gene_RRA$Name,up_gene_freq[,1]),2]
  # 计算下调基因的RRA排名
  down_gene_RRA<-aggregateRanks(glist=down_gene_list,N=length(unique(unlist(down_gene_list))))
  down_gene_freq<-as.data.frame(table(unlist(down_gene_list)))
  down_gene_RRA$Freq<-down_gene_freq[match(down_gene_RRA$Name,down_gene_freq[,1]),2]
  
  upg_sig_RRA=up_gene_RRA[up_gene_RRA$Score<0.05,1]
  downg_sig_RRA=down_gene_RRA[down_gene_RRA$Score<0.05,1]
  sig_RRA=upg_sig_RRA
}

RRA_pheatmap<-function(sig_RRA,gene_list_log2FC){
  gene_list4pheatmap<-NULL
  for( key in names(gene_list_log2FC)){

    item<-gene_list_log2FC[[key]]%>%filter(SYMBOL %in% sig_RRA)%>%dplyr::select(SYMBOL,log2FoldChange)
    if (is.null(gene_list4pheatmap)) {
      # 第一次迭代时初始化数据框
      gene_list4pheatmap <- item
      names(gene_list4pheatmap)[2]<-key
    }else{
        gene_list4pheatmap<-merge(gene_list4pheatmap, item, by = "SYMBOL", all = TRUE)
        names(gene_list4pheatmap)[ncol(gene_list4pheatmap)]<-key
      }
  }
  rownames(gene_list4pheatmap)<-gene_list4pheatmap$SYMBOL
  gene_list4pheatmap<-gene_list4pheatmap[,-1]
  order_data<-gene_list4pheatmap[match(sig_RRA,rownames(gene_list4pheatmap)),]
  my_palette<-colorRampPalette(c("blue","white","red"))(n=299)
  pheatmap(order_data,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color = my_palette,
           breaks = seq(-max(order_data,na.rm=TRUE),max(order_data,na.rm=TRUE),length.out=299)
           )
}

RRA_pheatmap(upg_sig_RRA,gene_list_log2FC)
RRA_pheatmap(downg_sig_RRA,gene_list_log2FC)
