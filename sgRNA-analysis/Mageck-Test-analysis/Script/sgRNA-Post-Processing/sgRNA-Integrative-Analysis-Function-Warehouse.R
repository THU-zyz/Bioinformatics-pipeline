##### sgRNA integrative analysis Function Warehouse#####

create_output_dir4sgRNA <- function(output_path){
  "
  Description:
  ------------
  根据给定的结果文件夹路径，检查其是否存在，如果不存在则创造该文件夹。
  并且在该文件夹下，创造sgrna和gene文件夹。
  
  Parameter:
  ---------
  output_path: 结果文件夹路径。
  
  Return:
  -------
  
  "
  # 检查路径是否存在
  if (!file.exists(output_path)) {
    # 如果不存在，创建路径
    dir.create(output_path, recursive = TRUE)
  }
  
  # 创建gene和sgrna文件夹
  gene_path <- file.path(output_path, "gene")
  sgrna_path <- file.path(output_path, "sgrna")
  
  # 检查并创建gene子路径
  if (!file.exists(gene_path)) {
    dir.create(gene_path, recursive = TRUE)
  }
  
  # 检查并创建sgrna子路径
  if (!file.exists(sgrna_path)) {
    dir.create(sgrna_path, recursive = TRUE)
  }
  # 返回创建的路径
  return(list(gene_path = gene_path, sgrna_path = sgrna_path))
}


ReadRRA_pvalue <- function(gene_summary, score = c("lfc", "rra")[1]){
  "
  Description:
  ------------
  ReadRRA_pvalue 用于从给定的基因表达数据摘要中提取并处理基因相关的统计信息。
  
  Parameters:
  -----------
  gene_summary: 基因相关的统计信息的文件路径，txt格式。
  score: 指定使用lfc(对数倍数变化)还是rra作为评分标准。
  
  Return:
  -------
  返回一个包含id,score,FDR等列的新数据框。
  "
  if(is.null(dim(gene_summary))){
    gene_summary = read.table(file = gene_summary, sep = "\t", header = TRUE, quote = "",
                              comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  }
  # 如果数据库已经存在id, score和FDR这三列，则直接使用这些列创建一个新的数据框dd并返回。
  if(all(c("id", "Score", "FDR")%in%colnames(gene_summary))){
    dd = as.data.frame(gene_summary[,c("id", "Score", "FDR")], stringsAsFactors = FALSE)
    dd$id = as.character(dd$id)
    return(dd)
  }
  # 选定指定列，并给他们赋予新的列名。
  gene_summary = gene_summary[, c(1, 3, 9, 8, 14, 5, 11,4,10,6,12)]
  colnames(gene_summary) = c("id", "negscore", "poscore", "neglfc", "poslfc", "negfdr", "posfdr","negPvalue","posPvalue","negRank","posRank")
  dd = gene_summary
  # 如果设置为lfc, 
  if("lfc" %in% tolower(score)){
    dd$LFC = dd$poslfc
    dd$FDR = dd$posfdr
    dd$Pvalue=dd$posPvalue
    # 对于负LFC大于正LFC的情况，使用负LFC对应值。
    dd$LFC[abs(dd$neglfc)>dd$poslfc] = dd$neglfc[abs(dd$neglfc)>dd$poslfc]
    dd$FDR[abs(dd$neglfc)>dd$poslfc] = dd$negfdr[abs(dd$neglfc)>dd$poslfc]
    dd$Pvalue[abs(dd$neglfc)>dd$poslfc] = dd$negPvalue[abs(dd$neglfc)>dd$poslfc]
    dd = dd[, c("id", "LFC", "FDR","Pvalue","negRank","posRank","negscore","poscore")]
  }else if("rra" %in% tolower(score)){
    # 计算每行的最大负对数评分，并根据此更新LFC和FDR列。
    # 如果负评分小于正评分，将LFC设置为负值。
    idx_neg = dd$negscore<dd$poscore
    dd$LFC = apply(-log10(dd[, 2:3]), 1, max);
    dd$LFC[idx_neg] = -dd$LFC[idx_neg]
    dd$FDR = dd$posfdr; dd$FDR[idx_neg] = dd$negfdr[idx_neg]
    dd = dd[, c("id", "LFC", "FDR")]
  }
  colnames(dd) = c("id", "Score", "FDR","Pvalue","negRank","posRank","negscore","poscore")
  dd$id = as.character(dd$id)
  return(dd)
}

sgRNA_analysis4gene<-function(gene_summary_path, output_path, gene_list, my_color_paltte, LFC_threshold, title_name){
  "
  Description:
  ------------
  sgRNA_analysis4gene用于对sgRNA的基因统计信息结果进行可视化分析，包括火山图分析。
  
  Parameters:
  -----------
  gene_summary_path: 基因相关的统计信息的文件路径，txt格式。
  output_path: 结果文件或图输出地址。
  gene_list: 基因列表，用于在图中展示这些基因。
  my_color_paltte: 自定义颜色盘。
  LFC_threshold: LFC阈值，当LFC_threshold>0时，选择正向筛选的结果；当LFC_threshold<0时，展示负向筛选的结果；
  title_name: 标题名称。
  
  Return:
  -------
  绘制并保存gene summary数据的Vocano plot和Rank plot。
  "
  gdata = ReadRRA_pvalue(gene_summary_path)
  # 筛选掉数据中以Intergenic或以Offtarget开头的id所在行，这些基因被认为是脱靶的或者是非有效的基因。
  gdata<-gdata%>%mutate(sig=ifelse(grepl("^Intergenic",id)|grepl("^Offtarget",id),"nontargeted","targeted"))
  gdata<-gdata%>%filter(sig!="nontargeted")
  # 计算FDR、Pvalue的-log10值。
  gdata$LogFDR = -log10(gdata$FDR)
  gdata$LogPvalue = -log10(gdata$Pvalue)
  
  ##### 绘制火山图 #####
  if(LFC_threshold <= 0){
    gdata<-gdata%>%mutate(color=ifelse(Pvalue>0.05,"gray",ifelse(Score > LFC_threshold,"gray","blue")))
  }else{
    gdata<-gdata%>%mutate(color=ifelse(Pvalue>0.05,"gray",ifelse(Score < LFC_threshold,"gray","red")))
  }
  
  p1<-gdata%>%ggplot(aes(x=Score,y=LogPvalue,color=color))+geom_point(size=2,alpha=0.8)+
    labs(x="Log2(Foldchange)",y="-Log10(P-Value)",title = title_name)+
    # 修改颜色
    scale_color_manual(values = c("gray" = "gray", "blue" = "#5599FF","red" = "#EE6363"))+
    theme(text = element_text(colour = "black", size = 18),
          plot.title = element_text(hjust = 0.5, size = 18),
          axis.text = element_text(colour = "gray10"),
          #axis.text.y = element_text(face = "italic"),
          panel.border = element_blank())+
    theme(axis.line = element_line(), 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          #panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"), 
          panel.background = element_blank())+
    theme(legend.position = "none")+
    xlim(-5,10)+
    geom_vline(xintercept =LFC_threshold,linetype = 3)+
    geom_hline(yintercept = -log10(0.05),linetype = 3)
  
  p1<-p1+annotate("text", x = -1, y = 4, label = sprintf("LFC = %s", LFC_threshold), size = 6) +
    annotate("text", x = 9, y = 1.5, label = "P-Value = 0.05", size = 6)
  
  p1<-p1+geom_text(data = subset(gdata, id %in% gene_list), aes(label = id), nudge_x = 0.1, nudge_y = 0.1, size = 4, color = "black")
  # 存储火山图结果文件
  ggsave(p1,file=paste0(output_path,'gene/',title_name," Volcano_plot.pdf"),width=11.69,height=8.27)
  print(paste0("Complete gene summary Vocano plot of ",title_name))
  ##### 绘制Rank图 #####
  # 当LFC_threshold<=0时绘制负向选择的Rank图，否则绘制正向选择的Rank图。
  if(LFC_threshold <= 0){
    gdata$Rank = gdata$negRank
    gdata<-gdata%>%mutate(lognegscore=-log10(negscore))
    gdata<-gdata%>%mutate(col=ifelse(!(id %in% gene_list),"gray",id))
    p2<-gdata%>%ggplot(aes(x=Rank,y=lognegscore,color=col))+geom_point(alpha=1,size=3)+labs(x="Gene Rank",y="MAGeCK RRA Score (-Log10)",title = paste0( title_name," Negative Rank"))
  }else{
    gdata$Rank = gdata$posRank
    gdata<-gdata%>%mutate(logposcore=-log10(poscore))
    gdata<-gdata%>%mutate(col=ifelse(!(id %in% gene_list),"gray",id))
    p2<-gdata%>%ggplot(aes(x=Rank,y=logposcore,color=col))+geom_point(alpha=1,size=3)+labs(x="Gene Rank",y="MAGeCK RRA Score (-Log10)",title = paste0( title_name," Postive Rank"))
  }
  
  p2<-p2+theme(text = element_text(colour = "black", size = 18),plot.title = element_text(hjust = 0.5, size = 18),axis.text = element_text(colour = "gray10"),panel.border = element_blank())
  p2<-p2+theme(axis.line = element_line(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.background = element_blank())
  # 修改颜色
  p2<-p2+scale_color_manual(values = c("gray" = "gray",my_color_paltte),breaks = my_color_paltte)
  
  
  p2<-p2+geom_text_repel(data = subset(gdata, id %in% gene_list),aes(label = id),size = 4,color = "black",box.padding = 0.5,nudge_x = 60)
  
  if(LFC_threshold <= 0){
    ggsave(p2,file=paste0(output_path,'gene/',title_name," Negative_Rank.pdf"),width=11.69,height=8.27)
    print(paste0("Complete gene summary Negative Rank plot of ",title_name))
  }else{
    ggsave(p2,file=paste0(output_path,'gene/',title_name," Postive_Rank.pdf"),width=11.69,height=8.27)
    print(paste0("Complete gene summary Postive Rank plot of ",title_name))
  }
  
}


sgRNA_analysis4sgrna<-function(sgrna_summary_path,gene_summary_path, geneset, output_path, title_name, binwidth = 0.1,interval = 0){
  "
  Description:
  -----------
  sgRNA_analysis4sgrna 用于可视化指定基因集中每一个基因的所有sgRNA在所有sgRNA中的LFC变化分布，从而帮助确定是哪一块的sgRNA起到关键的作用。
  
  Parameters:
  -----------
  sgrna_summary_path: sgrna summary文件所在路径。
  gene_summary_path: gene summary 文件所在路径。
  geneset: 指定可视化的基因集。
  output_path: 设置输出结果路径。
  title_name: 设置文件名称。
  binwidth,interval:画图参数，用于设置范围。
  Returns:
  --------
  "
  # 读取sgRNA_summary与gene_summary文件
  df = read.table(sgrna_summary_path, sep='\t', header=TRUE)
  df_gene = read.table(gene_summary_path, sep='\t', header=TRUE, row.names=1)
  # 选择df中的sgrna、Gene、LFC以及FDR列
  df <- df[,c('sgrna', 'Gene', 'LFC', 'FDR')]
  # 由于数据中包含一部分非记忆的sgRNA片段，在正常分析流程中这一部分片段将作为背景。
  # 匹配Gene列中以Intergenic开头的Gene,将其的sig列标记为nontargeted
  df_control<-df%>%mutate(sig=ifelse(grepl("^Intergenic",Gene),"nontargeted","targeted"))
  # 将sig=='nontargeted'的数据提取出来作为control.
  df_control <- df_control[df_control[,'sig']=='nontargeted',]
  df_control <- df_control[,c('sgrna', 'Gene', 'LFC', 'FDR')]
  
  # 筛选出指定的geneset
  subdf_geneset<-df[df[,'Gene']%in%geneset,]
  # 为df_control、subdf_geneset匹配上相应的颜色
  df_control$color <- "ctrl"
  subdf_geneset$color <- "Down"
  # 如果有多个subdf则用rbind将它们合并在一起。
  subdf<-subdf_geneset
  
  Gene <- unique(subdf$Gene)
  subdf$Gene <- factor(subdf$Gene, levels = geneset)
  subdf <- subdf[order(subdf$Gene), ]
  
  # 根据Gene的长度来扩增df_control的长度。
  df_control <- do.call("rbind", replicate(length(Gene), df_control, simplify = FALSE))
  
  # 计算画图参数，
  subdf$index <- rep(1:length(Gene), as.numeric(table(subdf$Gene)[Gene]))
  subdf$yend <- (binwidth + interval) * subdf$index - interval
  subdf$y <- (binwidth + interval) * (subdf$index - 1)
  
  df_control$index <- rep(1:length(Gene), rep(length(unique(df_control$sgrna)), length(Gene)))
  df_control$yend <- (binwidth + interval) * df_control$index - interval
  df_control$y <- (binwidth + interval) * (df_control$index - 1)
  
  # 对subdf的FDR进行整合，当lfc<0时,取neg.fdr的值，当lfc>0时取pos.fdr的值。
  subdf$all_FDR <- 0
  for (gene in Gene) {
    gene_lfc <- df_gene[gene,'neg.lfc']
    if (gene_lfc < 0) {
      gene_fdr <- df_gene[gene,'neg.fdr']
    } else {
      gene_fdr <- df_gene[gene,'pos.fdr']
    }
    
    subdf[subdf[,'Gene']==gene,'all_FDR'] <- gene_fdr
  }
  
  # 对subdf的p.value进行整合，当lfc<0时，取neg.p.value的值；
  # 当lfc>0时，取pos.p.value的值。
  subdf$all_pvalue <- 0
  for (gene in Gene) {
    gene_pvalue <- df_gene[gene,'neg.p.value']
    if (df_gene[gene,"neg.lfc"] < 0) {
      gene_pvalue <- df_gene[gene,'neg.p.value']
    } else {
      gene_pvalue <- df_gene[gene,'pos.p.value']
    }
    
    subdf[subdf[,'Gene']==gene,'all_pvalue'] <- gene_pvalue
  }
  
  a <- -Inf
  b <- Inf
  
  bindex <- as.vector(sapply(seq(1, max(subdf$index), 1), 
                             function(x) {
                               rep(x, 4)
                             }))
  bgcol <- data.frame(as.vector(bindex))
  bgcol$color <- c(rep("bg", length(bindex)))
  colnames(bgcol) <- c("id", "value")
  bgcol$x <- rep(c(a, b, b, a), max(subdf$index))
  bgcol$y <- as.vector(sapply(seq(1, max(subdf$index), 1), 
                              function(x) {
                                c((interval + binwidth) * (x - 1), (interval + binwidth) * 
                                    (x - 1), (interval + binwidth) * x - interval, 
                                  (interval + binwidth) * x - interval)
                              }))
  
  # 这里选择pvalue作为筛选条件
  label_list <- subdf[,c("Gene", "index", "all_pvalue")]
  label_list <- label_list[order(label_list$index,decreasing = F),]
  tmp_lab_dat <- as.data.frame(matrix(1:ncol(label_list),nrow = 1))
  colnames(tmp_lab_dat) <- colnames(label_list)
  tmp_lab_dat <- tmp_lab_dat[-1,]
  for (i in unique(label_list$index)) {
    ind_dat  <- label_list[label_list$index==i,]
    ind_dat <- ind_dat[1,]
    tmp_lab_dat <- rbind(tmp_lab_dat,ind_dat)  
  }
  label_list <- tmp_lab_dat
  # 为pvalue 贴上显著性标签，当pvalue>0.05为ns,当pvalue<=0.05为*,当pvalue<=0.01为**，当pvalue<0.001为***。             
  label_list$pvalue_label <- "a"
  for (i in 1:nrow(label_list)) {
    if (label_list[i,'all_pvalue']>0.05) {
      label_list[i,'pvalue_label'] <- "ns"
    }else if (label_list[i,'all_pvalue']<=0.05 & label_list[i,'all_pvalue']>0.01){
      label_list[i,'pvalue_label'] <- "*"
    }else if (label_list[i,'all_pvalue']<=0.01 & label_list[i,'all_pvalue']>0.001){
      label_list[i,'pvalue_label'] <- "**"
    }else{
      label_list[i,'pvalue_label'] <- "***"
    }
  }
  
  label_list <- label_list$pvalue_label
  # 为不同类型的基因设置不同颜色。
  cols <- c(Up = "#E41A1C", Down = "#377EB8", ctrl = "#4DAF4A" , tbg = 608, black = "black")
  
  # 绘制条形码图，其中y轴上，每一行代表一个基因，x轴为LFC的分布轴，每一个条带上的竖线代表这个基因对应的sgRNA在LFC上所在位置
  p <- ggplot()
  p <- p + geom_polygon(aes_string("x", "y", fill = "value", 
                                   group = "id"), color = "gray20", data = bgcol, size = 0.5)
  #p <- p + geom_segment(aes_string("LFC", "y", xend = "LFC", 
  #                                 yend = "yend", color = "color"), size = 0.5, alpha = 0.1, data = df_control)
  p <- p + geom_segment(aes_string("LFC", "y", xend = "LFC", yend = "yend", color = "color"), size = 0.8, data = subdf)
  p
  p <- p + scale_color_manual(values = cols)
  p <- p + scale_fill_manual(values = c(bg = 'white'))
  
  p <- p + scale_y_continuous(breaks = bgcol$y[seq(1, nrow(bgcol),4)] + binwidth/2,
                              labels = geneset,
                              expand = c(0, 0))
  
  p <- p + labs(x = 'LFC', y = NULL)
  # 设置坐标轴字体、颜色
  p <- p + theme(text = element_text(colour = "black", size = 18),
                 plot.title = element_text(hjust = 0.5, size = 18),
                 axis.text = element_text(colour = "gray10"),
                 axis.text.y = element_text(face = "italic"),
                 axis.ticks.y = element_blank(),)
  
  p <- p + theme(axis.line = element_blank(), panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"), 
                 panel.background = element_blank())
  # 设置是否设立图例
  p <- p + theme(legend.position = "none")
  # 设置辅助线
  p <- p + geom_vline(xintercept =c(0),linetype = 3)
  # p <- p + xlim(-5,5)
  # 设置显示范围。
  p <- p + xlim(-11,11)
  # 绘制sgRNA显著性标志
  panel_ann <- ggplot()+ 
    annotate("text",x = 0,y = bgcol$y[seq(1, nrow(bgcol),4)] + binwidth/2,label =label_list)  +
    theme_void() +
    theme(text = element_text(colour = "black", size = 18)) +
    xlab("") +ylab("") + scale_y_continuous(breaks = bgcol$y[seq(1, nrow(bgcol),4)] + binwidth/2,label=rep("",length(label_list)), expand = c(0, 0)) 
  
  df_noctrl<-df%>%mutate(sig=ifelse(grepl("^Intergenic",Gene)|grepl("^Offtarget",Gene),"nontargeted","targeted"))
  df_noctrl<-df_noctrl%>%filter(sig!="nontargeted")
  # 绘制sgRNA分布密度图
  density_plot <- ggplot(df_noctrl,aes(x = LFC)) +
    geom_density(color = "black",fill = "gray") +
    xlab("") + ylab("") + 
    # ggtitle('PANC02 E:T=1:1') +
    ggtitle('Negative regulators') +
    theme_bw() +
    theme(panel.border = element_rect(fill=NA, color="black", size=1, linetype="solid"),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    annotate("text", x = Inf, y = Inf, label = 'All sgRNAs', vjust = 1.8, hjust = 1.4)
  
  p <- p %>% insert_right(panel_ann,width = 0.1) %>% insert_top(density_plot,height = 0.15)
  ggsave(p,file=paste0(output_path,'sgrna/',title_name," Negative_regulators.pdf"),width=11.69,height=8.27)
}




