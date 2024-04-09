##### ssGSEA function repo #####
# read4ENSEMBL2SYMBOL读取count matrix并将其转换成以SYMBOL为行名的矩阵。
read4ENSEMBL2SYMBOL<-function(exp){
  exp<-rename(exp,ENSEMBL=gene_id)
  gene_symbol<-clusterProfiler::bitr(exp$ENSEMBL,fromType = "ENSEMBL",toType = "SYMBOL",OrgDb = "org.Hs.eg.db")
  exp<-exp%>%inner_join(gene_symbol,by="ENSEMBL")
  exp<-exp%>%select(SYMBOL,everything(),-ENSEMBL)
  # 使用dplyr 中的group_by和summarize_at函数对除了SYMBOL列外的所有列进行加和。
  df <- exp%>%
    group_by(SYMBOL) %>%
    summarize(across(everything(), sum))
  # 用SYMBOL作为行名并去除SYMBOL列。
  df<-as.data.frame(df)
  rownames(df)<-df$SYMBOL
  df<-df%>%select(-SYMBOL)
  # 去除表达谱中在所有样本中均未表达的基因。
  df<-df[rowSums(df[,-1])!=0,]
  # 将表达谱的data.frame 转化为matrix，否则后续的gsva分析会报错。
  df<-as.matrix(df)
  return(df)
}





