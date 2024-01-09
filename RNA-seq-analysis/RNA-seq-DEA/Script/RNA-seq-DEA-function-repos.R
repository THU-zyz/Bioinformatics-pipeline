##### RNA-seq differential expression analysis function repos #####

ENSEMBL2SYMBOL<-function(gene_matrix,OrgDb = "org.Hs.eg.db"){
  gene_SYMBOL<-bitr(gene_matrix$ENSEMBL, fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = "org.Hs.eg.db")
  df<-gene_matrix%>%inner_join(gene_SYMBOL,by="ENSEMBL")
  df<-df%>%dplyr::select(SYMBOL,everything(),-ENSEMBL)
  return(df)
}


volcano_plot<-function(res,output_name,plot_title,xlab,ylab){
  ###  Figure.1 Volcano plot illustrating genome-wide differential gene expression. ### 
  # Create a scatter plot;
  p12<-ggplot(data=res,aes(x=log2FoldChange,y=-log10(padj),color=sig))+geom_point(size=1)
  # Set point colors
  p12<-p12+scale_color_manual(values=c("red","gray","blue"),limits=c('up','none','down'))
  # Define axis titles
  p12<-p12+labs(x=xlab,y=ylab,title = plot_title)
  # Set background and font styles.
  p12<-p12+theme(plot.title = element_text(hjust=0.5,size=14),panel.grid = element_blank(),
                 panel.background = element_rect(color='black',fill='transparent'),
                 legend.key = element_rect(fill='transparent'),legend.title = element_blank(),
                 axis.title.x =element_text(size=12),axis.title.y=element_text(size=12)) 
  # Set threshold lines and axis range lines.
  p12<-p12+geom_vline(xintercept = c(-L2FCThread,L2FCThread),lty = 3, color = 'black')+
    geom_hline(yintercept = -log10(padjThread),lty=3,color='black')
  
  # Save the image as a PDF.
  pdf(paste0("Result/",output_name,"_volcano_plot.pdf"),width=6,height = 5)
  print(p12)
  dev.off()
}

Differential_expression_analysis<-function(input_files,sample_info,comparisons,output_name,human_metabolic_gene,plot_title,xlab,ylab,padjThread=0.05,readssumThread=1,L2FCThread=1.5){
  # read the input data.
  raw_count<-read.table(input_files,sep='\t',header=T)
  row.names(raw_count)<-raw_count[,1]
  raw_count <-raw_count[,-1]
  # Filter out genes with a total read count of less than 1 in all samples to obtain
  # the gene expression matrix.
  countdata<-raw_count[rowSums(raw_count)>=readssumThread,]
  # Based on the comparison and sample information, select all samples included in the two groups for differential analysis analysis.
  sample_info<-sample_info%>%filter(group %in% comparisons)
  countdata<-countdata[,sample_info$sample]
  ### 2.1) provide sample information
  condition_merge<-factor(sample_info$group)
  ### 2.2) Associate grouping information with sample names.
  colData<-data.frame(row.names=sample_info$sample,condition_merge)
  ### 3) DESeq analysis
  ### 3.1) Create a DESeq2Dataset object.
  dds<-DESeqDataSetFromMatrix(countdata,colData,design = ~condition_merge)
  ### 3.2) Normalization, involving multiple steps.
  dds2<-DESeq(dds)
  ### 4) Retrieve the results.
  res<-results(dds2)
  res<-res[order(res$padj),]
  res<-as.data.frame(res)
  # Generate a new column 'sig' based on 'padj' and 'log2FoldChange' value, assigning information on
  # gene upregulation or downregulation.
  res<-res %>%mutate(
    sig = case_when(
      padj < padjThread & log2FoldChange >= L2FCThread ~ "up",
      padj < padjThread & log2FoldChange <= -L2FCThread ~ "down",
      TRUE ~ "none"
    )
  )
  res$ENSEMBL<-rownames(res)
  res<-ENSEMBL2SYMBOL(res,OrgDb = "org.Hs.eg.db")
  # Save the differential expression results of all genes.
  write.csv(res,paste0("Result/",output_name,"-DESeq2-analysis-result.csv"))  
  # Filter out the differential expression genes.
  DE_res<-res%>%filter(sig!='none')
  write.csv(DE_res,paste0("Result/",output_name,"-DESeq2-analysis-DE-result.csv"))
  
  # Draw a volcano plot.
  volcano_plot(res,output_name,plot_title,xlab,ylab)
  
  # Draw a volcano plot for metabolism-related genes.
  geneInfo<-read.csv(human_metabolic_gene)  
  mr_res<-res%>%inner_join(geneInfo,by="SYMBOL")
  write.csv(mr_res,paste0("Result/",output_name,"-DESeq2-analysis-metabolic-related-gene-result.csv"))  
  # Filter out DE genes.
  DE_mr_res<-mr_res%>%filter(sig!='none')
  write.csv(DE_mr_res,paste0("Result/",output_name,"-DESeq2-analysis-metabolic-related-gene-DE-result.csv"))
  
  # Draw a volcano plot for metabolism-related genes
  volcano_plot(mr_res,paste0(output_name,"_metabolic_related_gene"),plot_title,xlab,ylab)
  
  results<-list(all=res,mr=mr_res)
  return(results)
}


Differential_expression_analysis<-function(input_files,sample_info,comparisons,refgroup,output_name,human_metabolic_gene,plot_title,xlab,ylab,padjThread=0.05,readssumThread=1,L2FCThread=1.5){
  # Read the input data
  raw_count<-read.table(input_files,sep='\t',header=T)
  row.names(raw_count)<-raw_count[,1]
  raw_count <-raw_count[,-1]
  # Filter out genes with a total read count of less than 1 in all samples to obtain the gene expression matrix.
  countdata<-raw_count[rowSums(raw_count)>=readssumThread,]
  # Based on the comparison and sample information, select all samples included in the two group for differential analysis.
  sample_info<-read.csv(sample_info_path)
  sample_info<-sample_info%>%filter(group %in% comparisons)
  countdata<-countdata[,sample_info$sample]
  ### 2.1) provide sample information
  condition_merge<-factor(sample_info$group)
  ### 2.2) Associate grouping information with sample names.
  colData<-data.frame(row.names=sample_info$sample,condition_merge)
  ### 3) DESeq analysis
  ### 3.1) Create a DESeq2Dataset object.
  dds<-DESeqDataSetFromMatrix(countdata,colData,design = ~condition_merge)
  dds$condition_merge<-relevel(dds$condition_merge,ref=refgroup)
  ### 3.2) Normalization, involving multiple steps.
  dds2<-DESeq(dds)
  ### 4) Retrieve the results.
  res<-results(dds2)
  res<-res[order(res$padj),]
  res<-as.data.frame(res)
  # Generate a new column 'sig' based on 'padj' and 'log2FoldChange' value, assigning information on gene upregulation or downregulation.
  res<-res %>%mutate(
    sig = case_when(
      padj < padjThread & log2FoldChange >= L2FCThread ~ "up",
      padj < padjThread & log2FoldChange <= -L2FCThread ~ "down",
      TRUE ~ "none"
    )
  )
  res$ENSEMBL<-rownames(res)
  res<-ENSEMBL2SYMBOL(res,OrgDb = "org.Hs.eg.db")
  # Save the differential expression results of all genes.
  write.csv(res,paste0("Result/",output_name,"-DESeq2-analysis-result.csv"))  
  # filter out DE gene
  DE_res<-res%>%filter(sig!='none')
  write.csv(DE_res,paste0("Result/",output_name,"-DESeq2-analysis-DE-result.csv"))
  
  # Draw a volcano plot.
  volcano_plot(res,output_name,plot_title,xlab,ylab)
  
  # Draw a volcano plot for metabolism-related genes 
  geneInfo<-read.csv(human_metabolic_gene)  
  mr_res<-res%>%inner_join(geneInfo,by="SYMBOL")
  write.csv(mr_res,paste0("Result/",output_name,"-DESeq2-analysis-metabolic-related-gene-result.csv"))  
  # filter out DE gene
  DE_mr_res<-mr_res%>%filter(sig!='none')
  write.csv(DE_mr_res,paste0("Result/",output_name,"-DESeq2-analysis-metabolic-related-gene-DE-result.csv"))
  
  # Draw a volcano plot for metabolism-related genes.
  volcano_plot(mr_res,paste0(output_name,"_metabolic_related_gene"),plot_title,xlab,ylab)
  
  results<-list(all=res,mr=mr_res)
  return(results)
}

#################### Customize the 'KEGGpathwayEnrich' function for pathway enrichment calculation. ####################
KEGGpathwayEnirch<-function(gene_set,fromType,toType,OrgDb="org.Hs.eg.db",organism="hsa",pvalueCutoff = 0.05,qvalueCutoff = 0.05){
  ##### Use 'bitr' for conversion based on the data type of the gene set. #####
  ##### fromType acceptes the gene type needed for conversion, while toType returns the gene type to be converted to, OrgDb accepts the reference database org.Hs.eg.db #####
  ##### The data commonly used for gene enrichment analysis is typically in the form of ENTREZID. #####
  gene_ENTREZID<-bitr(gene_set,fromType = fromType,toType = toType,OrgDb = OrgDb)
  kegg<-enrichKEGG(gene_ENTREZID$ENTREZID,organism = organism,
                   keyType = 'kegg',
                   pvalueCutoff = pvalueCutoff,
                   pAdjustMethod = 'BH',
                   qvalueCutoff = qvalueCutoff,
                   use_internal_data = FALSE)
  kegg<-setReadable(kegg,OrgDb = OrgDb,keyType = toType)
  kegg<-data.frame(kegg)
  ##### Calculate FoldEnrich based on GeneRatio and BgRatio #####
  kegg$FoldEnrich <- apply(kegg, 1,
                           function(x) as.numeric(unlist(strsplit(x[3], '[/]'))[1])/
                             as.numeric(unlist(strsplit(x[3], '[/]'))[2])*
                             as.numeric(unlist(strsplit(x[4], '[/]'))[2])/
                             as.numeric(unlist(strsplit(x[4], '[/]'))[1]))
  return(kegg)
}

# Input KEGG enrichment results, output KEGG_enrich_plot
KEGG_enrich_plot<-function(kegg_table,output_name){
  # Calculate the -log10FDR of KEGG_enrich_plot and the range of FoldEnrich。
  min_log10FDR<-floor(min(-log10(kegg_table$p.adjust),na.rm=T))
  max_log10FDR<-ceiling(max(-log10(kegg_table$p.adjust),na.rm = T))
  min_FoldEnrich<-floor(min(kegg_table$FoldEnrich,na.rm=T))
  max_FoldEnrich<-ceiling(max(kegg_table$FoldEnrich,na.rm=T))
  # Sort the data in ascending order based on FoldEnrich.
  kegg_table<-kegg_table[order(kegg_table$p.adjust,decreasing = T),]
  kegg_table$Description<-factor(kegg_table$Description,levels=unique(kegg_table$Description))
  # Draw a single-set KEGG enrichment diagram.
  p1<-ggplot(kegg_table,aes(x=-log10(p.adjust),y=Description))+
    geom_point(aes(fill = -log10(p.adjust),  size = FoldEnrich), shape = 21, color = "grey40")+
    scale_fill_viridis_c(direction = -1, end = 0.9, option = "C", limit = c(min_log10FDR, max_log10FDR)) + 
    scale_size(range = c(min_FoldEnrich, max_FoldEnrich)) + 
    theme_bw() +
    theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 2.5) +
    labs(x = "", y = "", title = "", size = "Fold enrichment", fill = "-Log10 (FDR)") +
    theme(panel.grid.major.x = element_line(size = 0), 
          panel.grid.minor.x = element_line(size = 0)) +
    theme(plot.title = element_text( size = 16, vjust = 6, hjust = 0.5, face = "plain"), 
          axis.title.x = element_text(size = 16, vjust = -2, hjust = 0.55, face = "plain"), 
          axis.title.y = element_text(size = 16, vjust = 2, face = "plain"),
          axis.text.x = element_text(size = 16, face = "plain", colour = "black"), 
          axis.text.y = element_text(size = 16, face = "plain", colour = "black"),
          legend.title = element_text(size = 12, face = "plain", colour = "black"),
          legend.text = element_text(size = 12), 
          legend.key.size = unit(0.4, "cm"), 
          legend.spacing = unit(0.4, "cm"))
  pdf(paste0("Result/",output_name,"-KEGG-Enrichmet-plot.pdf"),width=12,height=12)
  print(p1)
  dev.off()
}

# Input KEGG enrichment results, output KEGG pathways related to metabolic pathways
KEGG_enrich_metabolism_related_plot<-function(kegg_table,output_name,pathway_list_path){
  # Filter out the metabolic pathways among them
  pathway_list<-read.csv(pathway_list_path)
  Metabolism_pathway_list<-pathway_list%>%filter(Class_III=="Metabolism")
  kegg_table<-kegg_table%>%filter(ID %in% c(Metabolism_pathway_list$ID))
  # Calculate the -log10FDR of KEGG_enrich_plot and the range of FoldEnrich.
  min_log10FDR<-floor(min(-log10(kegg_table$p.adjust),na.rm=T))
  max_log10FDR<-ceiling(max(-log10(kegg_table$p.adjust),na.rm = T))
  min_FoldEnrich<-floor(min(kegg_table$FoldEnrich,na.rm=T))
  max_FoldEnrich<-ceiling(max(kegg_table$FoldEnrich,na.rm=T))
  # Sort the data in ascending order based on FoldEnrich
  kegg_table<-kegg_table[order(kegg_table$p.adjust,decreasing = T),]
  kegg_table$Description<-factor(kegg_table$Description,levels=unique(kegg_table$Description))
  # Draw a single-set KEGG Enrichment Diagram.
  p1<-ggplot(kegg_table,aes(x=-log10(p.adjust),y=Description))+
    geom_point(aes(fill = -log10(p.adjust),  size = FoldEnrich), shape = 21, color = "grey40")+
    scale_fill_viridis_c(direction = -1, end = 0.9, option = "C", limit = c(min_log10FDR, max_log10FDR)) + 
    scale_size(range = c(min_FoldEnrich, max_FoldEnrich)) + 
    theme_bw() +
    theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 2.5) +
    labs(x = "", y = "", title = "", size = "Fold enrichment", fill = "-Log10 (FDR)") +
    theme(panel.grid.major.x = element_line(size = 0), 
          panel.grid.minor.x = element_line(size = 0)) +
    theme(plot.title = element_text( size = 16, vjust = 6, hjust = 0.5, face = "plain"), 
          axis.title.x = element_text(size = 16, vjust = -2, hjust = 0.55, face = "plain"), 
          axis.title.y = element_text(size = 16, vjust = 2, face = "plain"),
          axis.text.x = element_text(size = 16, face = "plain", colour = "black"), 
          axis.text.y = element_text(size = 16, face = "plain", colour = "black"),
          legend.title = element_text(size = 12, face = "plain", colour = "black"),
          legend.text = element_text(size = 12), 
          legend.key.size = unit(0.4, "cm"), 
          legend.spacing = unit(0.4, "cm"))
  pdf(paste0("Result/",output_name,"-Metabolism-related-KEGG-Enrichmet-plot.pdf"),width=12,height=12)
  print(p1)
  dev.off()
}




multi_group_KEGG_enrich_plot<-function(mg_kegg_enrich,output_name){
  # Calculate the -log10FDR and the range of FoldEnrich for KEGG_enrich_plot.
  min_log10FDR<-floor(min(-log10(mg_kegg_enrich$p.adjust),na.rm=T))
  max_log10FDR<-ceiling(max(-log10(mg_kegg_enrich$p.adjust),na.rm = T))
  min_FoldEnrich<-floor(min(mg_kegg_enrich$FoldEnrich,na.rm=T))
  max_FoldEnrich<-ceiling(max(mg_kegg_enrich$FoldEnrich,na.rm=T))
  
  p2<-mg_kegg_enrich%>%ggplot(aes(x=cluster,y=reorder(Description,-p.adjust)))
  # Add a scatter plot with fill color
  p2<-p2+geom_point(aes(fill=-log10(p.adjust),size=FoldEnrich),shape=21,color="grey40")
  # Set the gradient direction for fill color, using the viridis color palette with the 'C' option
  p2<-p2+scale_fill_viridis_c(direction = -1, end = 0.9, option = "C",)
  # Set the theme to have a white background and remove panel grid lines
  # Adjust the plot margin and aspect ratio
  p2<-p2+theme_bw()+theme(panel.grid=element_blank())+
    theme(plot.margin = margin(1, 1, 1, 1, "cm"), aspect.ratio = 2.5)
  # Set labels for x-axis, y-axis, title, and legend with appropriate text sizes
  p2<-p2+labs(x = "", y = "", title = "", size = "Fold enrichment", fill = "-Log10 (FDR)")
  # Further customize the theme by removing major and minor grid lines on the x-axis
  # Customize font sizes and styles for title, axis titles, axis labels, and legend
  p2<-p2+theme(panel.grid.major.x = element_line(size = 0), 
               panel.grid.minor.x = element_line(size = 0))+
    theme(plot.title = element_text( size = 32, vjust = 6, hjust = 0.5, face = "plain"), 
          axis.title.x = element_text(size = 32, vjust = -2, hjust = 0.55, face = "plain"), 
          axis.title.y = element_text(size = 32, vjust = 2, face = "plain"),
          axis.text.x = element_text(size = 32, face = "plain", colour = "black"), 
          axis.text.y = element_text(size = 24, face = "plain", colour = "black"),
          legend.title = element_text(size =24, face = "plain", colour = "black"),
          legend.text = element_text(size = 24), 
          legend.key.size = unit(0.8, "cm"), 
          legend.spacing = unit(0.8, "cm"))
  p2<-p2+scale_fill_viridis_c(direction = -1, end = 0.9, option = "C", limit = c(min_log10FDR, max_log10FDR))
  p2<-p2+scale_size(range = c(4, 16))
  p2<-p2+scale_fill_gradientn(colors = c("#d2e032", "#31a881", "#3e5189"))
  pdf(paste0("Result/",output_name,"-KEGG-Enrichment-plot.pdf"),width=30,height=30)
  print(p2)
  dev.off()
}

multi_group_KEGGpathwayEnrich<-function(input_files,sample_info_path,
                                        comparisons_dir,refgroup,outputname_dir,plot_name,
                                        human_metabolic_gene,plot_title,xlab,ylab,
                                        padjThread=0.05,readssumThread=1,L2FCThread=1.5){
  all_kegg_enrich_list<-list()
  mr_kegg_enrich_list<-list()
  for(key in names(comparisons_dir)){
    comparisons <- comparisons_dir[[key]]
    output_name <- outputname_dir[[key]]
    allNmr<-Differential_expression_analysis(input_files,sample_info_path,
                                             comparisons,refgroup,output_name,human_metabolic_gene,plot_title,xlab,ylab,
                                             padjThread=0.05,readssumThread=1,L2FCThread=1.5)
    all_gene<-allNmr$all
    mr_gene<-allNmr$mr
    # Filter out differentially expressed genes
    all_DE_gene<-all_gene%>%filter(sig!="none")
    mr_DE_gene<-mr_gene%>%filter(sig!="none")
    # Perform KEGG enrichment for all genes
    all_kegg_enrich<-KEGGpathwayEnirch(gene_set=all_DE_gene$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db",organism="hsa",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
    # Perform KEGG enrichment for metabolism-related genes
    mr_kegg_enrich<-KEGGpathwayEnirch(gene_set=mr_DE_gene$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db",organism="hsa",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
    # Perform KEGG enrichment for all genes
    all_kegg_enrich<-KEGGpathwayEnirch(gene_set=all_DE_gene$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db",organism="hsa",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
    write.csv(all_kegg_enrich,paste0("Result/",output_name,"-KEGG-Enrichment-Result-All-gene.csv"),row.names = F)
    # Perform KEGG enrichment for metabolism-related genes
    mr_kegg_enrich<-KEGGpathwayEnirch(gene_set=mr_DE_gene$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db",organism="hsa",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
    write.csv(mr_kegg_enrich,paste0("Result/",output_name,"-KEGG-Enrichment-Result-mr-gene.csv"),row.names = F)
    # Draw a bubble chart containing all the enriched pathways
    KEGG_enrich_plot(all_kegg_enrich,paste0(output_name,"-all-"))
    # Draw a bubble chart for metabolic pathways
    KEGG_enrich_metabolism_related_plot(all_kegg_enrich,paste0(output_name,"-all-"),pathway_list_path)
    # Draw a bubble chart containing all the enriched pathways
    KEGG_enrich_plot(mr_kegg_enrich,paste0(output_name,"-mr-"))
    # Draw a bubble chart for metabolic pathways
    KEGG_enrich_metabolism_related_plot(mr_kegg_enrich,paste0(output_name,"-mr-"),pathway_list_path)
    
    all_kegg_enrich<-all_kegg_enrich%>%mutate(cluster=output_name)
    mr_kegg_enrich<-mr_kegg_enrich%>%mutate(cluster=output_name)
    # Store the results in a list
    all_kegg_enrich_list[[key]] <- all_kegg_enrich
    mr_kegg_enrich_list[[key]] <- mr_kegg_enrich
    print(paste0("KEGG pathway enrichment analysis for group ",output_name," has been completed."))
  }
  
  # Merge all databases after the loop ends
  all_kegg_enrich_combined <- do.call("rbind", all_kegg_enrich_list)
  write.csv(all_kegg_enrich_combined,paste0("Result/",output_name,"-KEGG-Enrichment-Result-All-gene-combined.csv"),row.names = F)
  mr_kegg_enrich_combined <- do.call("rbind", mr_kegg_enrich_list)
  write.csv(mr_kegg_enrich_combined,paste0("Result/",output_name,"-KEGG-Enrichment-Result-mr-gene-combined.csv"),row.names = F)
  # Draw bubble plots for multiple sets of KEGG enrichment
  multi_group_KEGG_enrich_plot(all_kegg_enrich_combined,paste0(plot_name,"-all-"))
  multi_group_KEGG_enrich_plot(mr_kegg_enrich_combined,paste0(plot_name,"-mr-"))
}