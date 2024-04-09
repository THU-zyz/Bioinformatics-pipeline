##### sgRNA Data Processing function warehouse #####
# The sgRNA Log2Readcount accept sgRNA read count data,
# converts it into a log2 format, and ultimately counts the number of occurrences at different values.

sgRNA_Log2Readcount<-function(input_files,output_dir,sample_info_path,plot_name){
  # 读取样本信息文件
  sample_info<-read.xlsx(sample_info_path)
  plasmid.count<-read.table(input_files,sep='\t',header = T)
  # Perform log2 transformation on the data from the third column to the last column.
  plasmid.count[,3:ncol(plasmid.count)]<-log2(plasmid.count[,3:ncol(plasmid.count)]+1)
  # Round down the data to two decimal places.
  plasmid.count[,3:ncol(plasmid.count)]<-floor(plasmid.count[,3:ncol(plasmid.count)]*100)/100
  for(sample_id in sample_info$Sample_ID){
    data<-plasmid.count%>%dplyr::select(all_of(sample_id))
    # sample_id是一个包含列名的变量，用`!!`和`sym`将字符串转换为符号
    # 确保`group_by`根据变量的值而不是变量名本身来分组。
    data<-data%>%group_by(!!sym(sample_id))%>%count()
    # plot sgRNA Log2Readcount.
    plot<-data%>%ggplot(aes(x=!!sym(sample_id),y=n))+
      geom_area(colour="black",fill="blue",alpha = 0.2)+
      scale_x_continuous(expand=c(0,0))+
      scale_y_continuous(expand=c(0,0))+
      theme_bw()+
      theme(#legend.position = "none",
        #axis.text.x=element_text(angle=30,hjust=1),
        plot.title = element_text(hjust=0.5),
        text = element_text(size=10,face="bold"),
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        strip.text.y = element_blank(), 
        #strip.text = element_text(size=10,face="bold"),
        strip.background = element_blank(),
        axis.line = element_line(color = 'black'))+
      labs(x="Log2 Readcount",
           y="number of sgRNAs",fill="carb")+
      ggtitle(plot_name)
    # save the plot
    ggsave(plot,file=paste0(output_dir,sample_id,"-log2Readcount.pdf"),width = 8,height = 6)
  }
}


