## signautre score analysis function repos ##
library(tidyverse)
library(openxlsx)
library(ggplot2)
library(ggpubr)
library(clusterProfiler)

## input a gene expression matrix, where the column name of the column containing gene ENSEMBL IDs is labeled as "ENSEMBL".
## return a gene expression matrix that includes the transformed 'SYMBOL' column information.
ENSEMBL2SYMBOL<-function(gene_matrix,OrgDb = "org.Hs.eg.db"){
  gene_SYMBOL<-bitr(gene_matrix$ENSEMBL, fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = "org.Hs.eg.db")
  df<-gene_matrix%>%inner_join(gene_SYMBOL,by="ENSEMBL")
  df<-df%>%dplyr::select(SYMBOL,everything(),-ENSEMBL)
  return(df)
}

## gene_matrix_preprocessing function input gene matrix, gene list, sample information，gene effective length.
## gene matrix is a gene expression matrix, the first column is ENSEMBL IDs name ENSEMBL and other column are sample names.
## geneList records Gene SYMBOL information, with the column named 'SYMBOL'.Specifically, sometimes geneList my also have a `Weight` column recording the trend of gene changes.
## sample_info records sample name and group information.
## gene_efflen record effective length for every gene from hg38 genome. 
## return a gene expression that selected by geneList and normalized by z-score.
gene_matrix_preprocessing<-function(gene_matrix,geneList,sample_info,gene_efflen){
  # 1) convert raw counts to tpm
  gene_matrix<-count2tpm(gene_matrix,gene_efflen)
  # 2) convert ENSEMBL ID to SYMBOL 
  gene_matrix<-ENSEMBL2SYMBOL(gene_matrix,"org.Hs.eg.db")
  # 3) merge rows with the same SYMBOL name,summing their TPM values
  unique_gene_matrix<-gene_matrix%>%group_by(SYMBOL)%>%summarise(across(everything(),sum))
  # 4) select  columns for analysis based on sample information
  gene_expression<-unique_gene_matrix[,intersect(names(unique_gene_matrix),c("SYMBOL",sample_info$sample)),drop=FALSE]
  # 5) perform z-score normalization on each row.
  gene_expression[,-1]<-t(scale(t(gene_expression[,-1])))
  # 6) select relevant genes based on a gene list.
  # if the number of columns in geneList exceeds 1, it implies the presence of a "Weight" column, which requires additional processing.
  gene_expression<-gene_expression%>%inner_join(geneList,by="SYMBOL")
  if(length(geneList)>1){
    gene_expression%>%dplyr::select(SYMBOL,Weight,everything())
    # if there are more than 1 columns in geneList, multiply their expression values by the corresponding weights.
    gene_expression[,-c(1,2)]<-gene_expression[,-c(1,2)]*gene_expression$Weight
    gene_expression<-gene_expression%>%dplyr::select(-Weight)
  }
  return(gene_expression)
}


# signature_score accept gene_matrix and sample_info as input,
# return a sign_score dataframe with an average score of all gene expression(z-score) in gene list.
signature_score<-function(gene_matrix,sample_info){
  
  gene_matrix<-t(gene_matrix)
  colnames(gene_matrix)<-gene_matrix[1,]
  gene_matrix<-as.data.frame(gene_matrix[-1,])
  gene_matrix<-gene_matrix%>%mutate_all(as.numeric)
  gene_matrix$Means<-rowMeans(gene_matrix)
  gene_matrix$sample<-rownames(gene_matrix)
  
  gene_matrix<-gene_matrix%>%inner_join(sample_info,by="sample")
  sign_score<-gene_matrix%>%dplyr::select(sample,group,sign_score=Means)
  return(sign_score)
}  

# signature_score_plot function is used to generate a scatter plot with error bars based on sign_score and compare significance.
# signature_score_plot function requires inputs sign_score and my_palette.
# By default, it compares the significance between multiple groups and a reference background (the average of all groups): comparisons='.all.'.
# Alternatively, custom comparisons can be selected to compare differences between specified groups.
# The Welch two-sample t-test is used by default, with **p < 0.01 and ***p < 0.001 as labels.
# The function ultimately returns a well-drawn plot.


signature_score_plot<-function(sign_score,title,my_palette,comparisons=".all.",ylab="Signature Score",Test_Methods="t.test",hide_ns=T,sign_label="p.signif"){
  # Compute and calculate the mean and standard deviation for each group in sign_score, which will be used to generate error bars in the plot.
  summary_data<-sign_score%>%group_by(group)%>%summarize(mean_score=mean(sign_score),sd_score=sd(sign_score),se_score=sd(sign_score)/sqrt(length(sign_score)))
  # 1) Generate scatter plots, error bars, and crossbars.
  p1<-sign_score%>%ggplot(aes(x=group,y=sign_score,color=group))+
    geom_point(aes(color =group),size=1,position = position_jitterdodge(dodge.width=0.1,jitter.width = 0.2))+
    geom_errorbar(data=summary_data,
                  aes(x=group,y=mean_score,ymin=mean_score-se_score,ymax=mean_score+se_score,color=group),
                  position=position_dodge(width=0.4),
                  width=0.2)+
    geom_crossbar(data=summary_data,aes(x=group,y=mean_score,ymin=mean_score,ymax=mean_score),width=0.3)
  # 2) Plot significance markers.
  # The t.test() function in R defaults to assuming unequal variances and uses the Welch method to correct degrees of freedom.
  if(".all." %in% comparisons){
    # The reference group is set to ".all.", representing the average of all groups. When comparing based on patient grouping, there are numerous combinations for inter-group comparisons.
    # Therefore, it is possible to compare each group from multiple groups with the overall average, observing whether the signature score is overexpressed or underexpressed in different groups.
    compare_means(sign_score~group,data=sign_score,ref.group = comparisons,method = Test_Methods)
    
    p2<-p1+stat_compare_means(
      ref.group = comparisons,
      method = Test_Methods,
      hide.ns =hide_ns,
      label=sign_label,
      size=10
    )
  } else{
    # Alternatively, perform differential testing between any two groups, and use horizontal lines to indicate the tested pairs.
    p2<-p1+stat_compare_means(
      comparisons = comparisons,
      method = Test_Methods,
      hide.ns =hide_ns,
      label=sign_label,
      size=10
    )
  }
  
  # 3) Configure plotting parameters such as font, color, size, and other visual aspects.
  p3<-p2+labs(y=ylab)+theme_bw()+ theme(
    legend.position="none",
    panel.grid=element_blank(),
    panel.border=element_blank(),
    # axis.line = element_line(size=1, colour = "black"),
    strip.background = element_blank(),
    strip.text = element_text(color="black",family = "Arial",size = 26),
    # axis.line = element_line(color="black"),
    axis.line.x = element_blank(),
    axis.line.y = element_line(color="black"),
    legend.title = element_text( color="black",family = "Arial", size=22),
    legend.text= element_text( color="black",family = "Arial", size=22),
    plot.title = element_text(hjust = 0.5,size=30,face="bold"),
    axis.text.x = element_text(face="bold",color="black",family = "Arial",size = 18),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face="bold",  color="black", size=15),
    axis.title.x = element_blank(),
    axis.title.y = element_text(face="bold",color="black", size=24))+
    scale_color_manual(values=my_palette)+ggtitle(title)
  
  return(p3)
}
