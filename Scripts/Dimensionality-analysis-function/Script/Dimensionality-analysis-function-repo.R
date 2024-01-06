#####  t-SNE analysis function repo #####
# t-dsitributed Stochastic Neighbor Embedding. 
library(tidyverse)
library(Rtsne)
library(openxlsx)
library(tidyverse)

# DataLoader is used for reading data and processing it into a format that is convenient for subsequent analysis. 
# The input data needs to have metabolite names as row names and sample names as column names, and the sample names need to follow the format 'group-sample number'.
# The final output will be a DataFrame, with the first column being the group, row names as sample names, and column names as metabolite names.
DataLoader<-function(files,sheetnames){
  Data<-read.xlsx(files,sheet=sheetnames)
  Data<-as.data.frame(t(Data))
  colnames(Data)<-Data[1,]
  Data<-Data[-1,]
  Data$Group<-rownames(Data)
  Data<-Data%>%mutate(Class=str_extract(Group,"^[^-]+"))
  Data<-Data%>%dplyr::select(-'Group')
  Data<-Data[c("Class",setdiff(colnames(Data),"Class"))]
  Data[, 2:ncol(Data)] <- lapply(Data[, 2:ncol(Data)], function(x) as.numeric(as.character(x)))
  
  return(Data)
}

# plot_tsne creates a t-SNE plot based on the results from tsne.result.
# tsne.results, output files name and color palette are required to be provided.
plot_tsne<-function(tsne.result,output_name,my_palatte){
  
  p1<-ggplot(tsne.result,aes(x=V1,y=V2,fill=group))+
    geom_point(size = 6, shape = 21, color = "black") +
    
    scale_fill_manual(values = my_palatte) +
    stat_ellipse(level = 0.95, show.legend = T, linetype = 2, color = "grey20")+
    labs(x = "tSNE-1", y = "tSNE-2", fill = "Group") + 
    theme_bw() +
    theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 1) +
    theme(axis.title.x = element_text(size = 16, hjust = 0.5, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
          axis.text.x = element_text(size = 16, face = "bold", colour = "black")) +
    theme(legend.text = element_text(size = 16), 
          legend.title = element_text(size = 16,face = "bold"), 
          legend.spacing.y = unit(0.2, 'cm'), 
          legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'))
  
  ggsave(p1,device = cairo_pdf,
         path="Result/",
         filename = paste0(output_name,".pdf"),
         width = 8,height = 10,
         # measurement units
         units ="in")
  return(p1)
}

# Rtsne() by default performs duplicate data detection during execution, and will throw an error if duplicate samples are not removed in advance.
# The results of t-SNE are random, so it's necessary to set a fixed random seed in advance to ensure the reproducibility of the results.
tSNE_analysis<-function(files,sheetnames,output_name,my_palatte,perplexity=5,normalized=T,seed=10086){
  Data<-DataLoader(files,sheetnames)
  colnumber=length(colnames(Data))
  # ensure numeric type.
  df.PCA<-data.frame(Data[,2:colnumber],check.names = F)
  # Set random seed
  set.seed(seed)
  tsne_out<-Rtsne(
    df.PCA,
    normalize =normalized,
    # set the number of dimensions after dimensionality reduction
    dims =2,
    # Whether to perform PCA analysis on the input raw data and then use top N PC obtained from PCA for subsequent analysis.
    # For data whith higher dimensions,PCA can be used to reduce dimensionality and improve operational efficiency.
    pca=F,
    # The 'perspective' generally ranges between 5-50, estimating the number of neighbors for each point.
    # 3 * perplexity < nrow(X) - 1, 
    perplexity=perplexity,
    # Theta represents the trade-off between running speed and accuracy. The default is 0.5; 0 represents exact t-SNE.
    theta =0.01,
    # The max iteration.
    max_iter = 1000
  )
  
  # Extract plotting information from the t-SNE results.
  tsne.result<-as.data.frame(tsne_out$Y,row.names = rownames(Data))
  tsne.result$group<-Data$Class
  tsne.result$name<-rownames(tsne.result)
  write.csv(tsne.result,paste0("Result/",output_name,".csv"),row.names = F)
  plot_tsne(tsne.result,output_name,my_palatte)
}

# Run PCA analysis, save the relevant results and create a PCA plot
pca_analysis<-function(files,sheetnames,output_name,my_palatte){
  Data<-DataLoader(files,sheetnames)
  colnumber<-length(colnames(Data))
  # z-score normalization.
  # The prcomp function includes a normalization step, so there is no need for additional normalization.
  pca_result<-prcomp(Data[,2:colnumber],scale.=TRUE,center=TRUE)
  # Extract the principal components.
  pca_data<-as.data.frame(pca_result$x)
  pca_data$group<-Data$Class
  write.csv(pca_data,paste0("Result/",output_name,".csv"),row.names = F)
  variance_explained<-summary(pca_result)$importance[2,1:2]*100
  # create a PCA plot
  p1<-ggplot(pca_data,aes(x=PC1,y=PC2,fill=group))+
    geom_point(size = 6, shape = 21, color = "black")+
    scale_fill_manual(values = my_palatte) +
    xlab(paste0("PC1 (",round(variance_explained[1],2),"%)"))+
    ylab(paste0("PC2 (",round(variance_explained[2],2),"%)"))+
    labs(fill = "Group")+ggtitle("")+
    stat_ellipse(level = 0.95, show.legend = T, linetype = 2, color = "grey20")+
    theme_bw() +
    theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 1) +
    theme(axis.title.x = element_text(size = 16, hjust = 0.5, face = "bold"),
          axis.title.y = element_text(size = 16, face = "bold"), 
          axis.text.y = element_text(size = 16, face = "bold", colour = "black"),
          axis.text.x = element_text(size = 16, face = "bold", colour = "black")) +
    theme(legend.text = element_text(size = 16), 
          legend.title = element_text(size = 16,face = "bold"), 
          legend.spacing.y = unit(0.2, 'cm'), 
          legend.key.height = unit(0.3, 'cm'), legend.key.width = unit(0.3, 'cm'))
  ggsave(p1,device = cairo_pdf,
         path="Result/",
         filename = paste0(output_name,".pdf"),
         width = 8,height = 10,
         # measurement units
         units ="in")
  return(p1)
}

# Run PCA analysis, save the relevant results and create a PCA plot
plsda_analysis<-function(files,sheetnames,output_name){
  Data<-DataLoader(files,sheetnames)
  colnumber<-length(colnames(Data))
  # z-score normalization.
  # The plsda function includes a normalization step, so there is no need for additional normalization.
  plsda_result<-plsda(Data[,-1],Data$Class,ncomp=2)
  pdf(paste0("Result/",output_name,".pdf"),width=10,height=8)
  plotIndiv(plsda_result,comp=c(1,2),
            group=Data$Class,style='ggplot2',ellipse = T,
            size.xlabel = 20,size.ylabel = 20,size.axis = 20,pch=16,cex=5,
            ind.names = FALSE,lengend=TRUE,title='PLS-DA'
  )
  dev.off()
  
}

