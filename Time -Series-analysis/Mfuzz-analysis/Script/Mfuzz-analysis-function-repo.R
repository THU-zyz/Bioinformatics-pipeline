##### Mfuzz function repo #####

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

# Dataloader function 
MDataLoader<-function(files,sheetnames){
  # Data preprocessing
  Data<-read.xlsx(files,sheet=sheetnames)
  Data<-as.data.frame(t(Data))
  colnames(Data)<-Data[1,]
  Data<-Data[-1,]
  Data$Group<-rownames(Data)
  # The sample names in the original file need to be standard.
  # The '-' character separates group information and group signals.
  # Extract the information before the '-'.
  Data<-Data%>%mutate(Class=str_extract(Group,"^[^-]+"))
  Data<-Data%>%dplyr::select(-'Group')
  Data<-Data[c("Class",setdiff(colnames(Data),"Class"))]
  # Set the first column, named 'Class'
  Data[,1]<-as.factor(Data[,1])
  Class=levels(Data[,1])
  classnumber=Class
  result<-list(Data=Data,Class=Class)
  return(result)
}

TDataLoader<-function(expression_path,sheetnames,sample_info_path,gene_efflen_path){
  # read gene expression matrix
  expression<-read.table(expression_path,header = TRUE)
  expression<-expression%>%dplyr::select(ENSEMBL=gene_id,everything())
  # read sample information 
  sample_info<-read.csv(sample_info_path,header=T,)
  # read gene effective length
  gene_efflen<-read.csv(gene_efflen_path)
  gene_matrix<-count2tpm(expression,gene_efflen)
  gene_matrix<-ENSEMBL2SYMBOL(gene_matrix,"org.Hs.eg.db")
  # merge rows with the same SYMBOL name,summing their TPM values
  unique_gene_matrix<-gene_matrix%>%group_by(SYMBOL)%>%summarise(across(everything(),sum))
  # select samples based on sample information.
  gene_expression<-unique_gene_matrix[,intersect(names(unique_gene_matrix),c("SYMBOL",sample_info$sample)),drop=FALSE]
  gene_expression_t<-as.data.frame(t(gene_expression))
  colnames(gene_expression_t)<-gene_expression_t[1,]
  gene_expression_t<-gene_expression_t[-1,]
  gene_expression_t$sample<-rownames(gene_expression_t)
  gene_expression_t<-gene_expression_t%>%inner_join(sample_info,by="sample")
  Data<-gene_expression_t%>%dplyr::select(Class=group,everything(),-sample)
  Data[,1]<-as.factor(Data[,1])
  Data[, 2:ncol(Data)] <- lapply(Data[, 2:ncol(Data)], function(x) as.numeric(as.character(x)))
  col_sums<-apply(Data[,-1], 2, sum)
  selected_cols<-which(col_sums>=1)
  Data<-Data[,c(1,selected_cols+1)]
  Class=levels(Data[,1])
  classnumber=Class
  result<-list(Data=Data,Class=Class)
  return(result)
}

KW_test<-function(paras,files,sheetnames,output_name){
  if(paras[["type"]]=="Metabolomics"){
    type=paras[["type"]]
    DataNclassnumber<-MDataLoader(files,sheetnames)
  }else if(paras[["type"]]=="Transcriptomics"){
    type=paras[["type"]]
    sample_info_path<-paras[["sample_info_path"]]     
    gene_efflen_path<-paras[["gene_efflen_path"]]
    DataNclassnumber<-TDataLoader(files,sheetnames,sample_info_path,gene_efflen_path)
  }
  
  df.test<-DataNclassnumber$Data
  # Ensure that the values within are numeric.
  df.test[, 2:ncol(df.test)] <- lapply(df.test[, 2:ncol(df.test)], function(x) as.numeric(as.character(x)))
  # Class number.
  Class<-DataNclassnumber$Class
  classnumber<-length(Class)
  # setting the empty dataframe to store the mediated files
  DunnRes.df <- data.frame() 
  result.dunn.pvalue <- data.frame()
  result.dunn.FDR<-data.frame()
  KW.Dunn <-data.frame()
  
  for (i in c(2:(length(df.test)))) {
    KruWRes<-kruskal.test(df.test[[i]]~Class,data=df.test)
    DunnRes.df <- rbind(DunnRes.df, unlist(KruWRes[3]))
    DunnRes <-dunn.test(df.test[[i]],df.test$Class,
                        altp = TRUE, method="bh")
    result.dunn.pvalue <- rbind(result.dunn.pvalue, unlist(DunnRes[3]))
    result.dunn.FDR <- rbind(result.dunn.FDR, unlist(DunnRes[4]))
    KW.Dunn <- cbind(DunnRes.df,result.dunn.pvalue,result.dunn.FDR)
    
  }
  # Name the row and column headers of Kruskal-Wallis test results.
  colnames(KW.Dunn) <- c("p-value: all",paste0("p-value: ", unlist(DunnRes[5])),paste0("FDR: ", unlist(DunnRes[5])))
  rownames(KW.Dunn) <- colnames(df.test)[2:(length(df.test))]
  # Calculate the Benjamini-Hochberg correction value for the Kruskal-Wallis p-value.
  KW.Dunn$`FDR: all`<-p.adjust(KW.Dunn$`p-value: all`,method="BH")
  KW.Dunn<-KW.Dunn%>%dplyr::select(`p-value: all`,`FDR: all`,everything())
  
  ##### Among multiple groups, perform pairwise comparisons to calculate the Fold Change(FC) #####
  # Set the names for multiple group comparisons.
  comparison <- combn(classnumber, 2)
  for(i in 1:length(Class)){
    comparison[comparison==i]<-Class[i]
  }
  df.test.fc <- data.frame(df.test, check.names = F)
  df.test.fc[, 2:ncol(df.test.fc)] <- lapply(df.test.fc[, 2:ncol(df.test.fc)], function(x) as.numeric(as.character(x)))
  df.test.fc$Class <- df.test$Class
  df.test.fc$Class <- as.factor(df.test.fc$Class)
  ##### calculate foldchange #####
  for (i in 1:(length(comparison)/2)){
    FCname <- paste0("FC","-", comparison[,i][1], "/", comparison[,i][2])
    KW.Dunn[[FCname]] <- apply(df.test.fc[, 2:(length(df.test))], 2, 
                               function(x) 
                                 mean(x[which(df.test.fc$Class == comparison[,i][1])])/
                                 mean(x[which(df.test.fc$Class == comparison[,i][2])]))
    LOG2FCname <- paste0("LOG2FC","-", comparison[,i][1], "/", comparison[,i][2])
    KW.Dunn[[LOG2FCname]] <- log2(KW.Dunn[[FCname]])
  }
  # order by p-value.
  KW.Dunn<-KW.Dunn[order(KW.Dunn$`p-value: all`),]
  write.csv(KW.Dunn,paste0("Result/",type,"/",output_name,'-KW-test-result.csv'))
  result<-list(KW.Dunn=KW.Dunn,df.test=df.test)
  return(result)
}

ENSEMBL2SYMBOL<-function(gene_matrix,OrgDb = "org.Hs.eg.db"){
  gene_SYMBOL<-bitr(gene_matrix$ENSEMBL, fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = "org.Hs.eg.db")
  df<-gene_matrix%>%inner_join(gene_SYMBOL,by="ENSEMBL")
  df<-df%>%dplyr::select(SYMBOL,everything(),-ENSEMBL)
  return(df)
}

Mfuzz_KW<-function(paras,files,sheetnames,output_name,pvalue_threshold,labels,Mfuzz_number){
  KW.DunnNdf.test<-KW_test(paras,files,sheetnames,output_name)
  type=paras[["type"]]
  KW.Dunn<-KW.DunnNdf.test$KW.Dunn
  df.test<-KW.DunnNdf.test$df.test
  # selected feature for clustering based on p-value threshold
  sign_KW<-KW.Dunn%>%dplyr::filter(`p-value: all`<pvalue_threshold)
  Mfuzz_data<-df.test[,c("Class",rownames(sign_KW))]  
  Mfuzz_data<-Mfuzz_data%>%group_by(Class)%>%summarise(across(everything(), mean, na.rm = TRUE))
  Mfuzz_data<-t(Mfuzz_data)
  colnames(Mfuzz_data)<-Mfuzz_data[1,]
  Mfuzz_data<-as.data.frame(Mfuzz_data[-1,])
  # Mfuzz data preprocessing
  Mfuzz_data_matrix<-as.matrix(sapply(Mfuzz_data, as.numeric))
  rownames(Mfuzz_data_matrix)<-rownames(Mfuzz_data)
  num<-ncol(Mfuzz_data_matrix)
  # construct an object 
  Mfuzz_class <- new('ExpressionSet',exprs = Mfuzz_data_matrix)
  # remove genes with too little variation between samples based on standard deviation.
  Mfuzz_class <- filter.std(Mfuzz_class, min.std = 0)
  # standardise
  Mfuzz_class <- standardise(Mfuzz_class)
  # The number of Mfuzz cluster 
  # Fuzzy c-means clustering requires manually defining the number of clusters.
  n <-Mfuzz_number
  # Evaluate the optimal m value to prevent clustering of random data.
  m <- mestimate(Mfuzz_class)
  # clustering
  set.seed(123)
  cl <- mfuzz(Mfuzz_class, c = n, m = m)
  # output the result 
  dat.cluster <- as.data.frame(cl$cluster)
  dat.cluster$metabolite <- row.names(dat.cluster)
  dat.cluster <- dat.cluster[order(dat.cluster$`cl$cluster`),]
  colnames(dat.cluster) <- c("group", "metabolite")
  # store the Mfuzz result files.
  write.csv(dat.cluster,paste0("Result/",type,"/",output_name,"-cluster.csv"),row.names = F)
  # save mfuzz plot.
  pdf(paste0("Result/",type,"/",output_name,"-plot.pdf"),width=8,height = 6)
  mfuzz.plot2(Mfuzz_class, cl = cl, mfrow = c(2, 2), ax.col = "black", cex = 3, x11 = F,
              xlab = "Time",
              ylab = "Relative abundance",
              cex.main = 2,
              #lwd=2,
              centre = T, centre.col = "grey20", centre.lwd = 3,
              time.labels = labels,
              colo = "fancy",
              #colo = colorRampPalette(c("#e11b47", "#ce1020"))(100),
              min.mem = 0)
  dev.off()
  
  result_data<-as.data.frame(t(as.data.frame(Mfuzz_class)))
  result_data$metabolite<-rownames(Mfuzz_data)
  rownames(result_data)<-result_data$metabolite
  result_data<-result_data%>%inner_join(dat.cluster,by="metabolite")
  
  write.xlsx(result_data, paste0("Result/",type,"/",output_name,"-result.xlsx"),sheetName = "4cluster", rowNames = F, colNames = T, append = T)
  return(result_data)
}


mfuzz.plot2 <- function(eset,cl,mfrow=c(1,1),colo,min.mem = 0,time.labels,time.points,ylim.set=c(0,0),
                        xlab="Time",ylab="Expression changes",x11=TRUE,
                        ax.col="black",bg = "white",col.axis="black",col.lab="black",                   
                        col.main="black",col.sub="black",col="black",centre=FALSE,
                        centre.col="black",centre.lwd=2,
                        Xwidth=5,Xheight=5,single=FALSE,...){
  # function for plotting the clusters 
  clusterindex <- cl[[3]]
  memship <- cl[[4]]
  memship[memship < min.mem] <- -1 
  colorindex <- integer(dim(exprs(eset))[[1]])
  
  
  if (missing(colo)){
    colo <- c("#FF0000","#FF1800",  "#FF3000", "#FF4800", "#FF6000", "#FF7800", "#FF8F00",
              "#FFA700", "#FFBF00", "#FFD700", "#FFEF00", "#F7FF00", "#DFFF00", "#C7FF00",
              "#AFFF00", "#97FF00", "#80FF00", "#68FF00", "#50FF00", "#38FF00", "#20FF00",
              "#08FF00", "#00FF10", "#00FF28", "#00FF40", "#00FF58", "#00FF70", "#00FF87",
              "#00FF9F", "#00FFB7", "#00FFCF", "#00FFE7", "#00FFFF", "#00E7FF", "#00CFFF",
              "#00B7FF", "#009FFF", "#0087FF", "#0070FF", "#0058FF", "#0040FF", "#0028FF",
              "#0010FF", "#0800FF", "#2000FF", "#3800FF", "#5000FF", "#6800FF", "#8000FF",
              "#9700FF", "#AF00FF", "#C700FF", "#DF00FF", "#F700FF", "#FF00EF", "#FF00D7",
              "#FF00BF", "#FF00A7", "#FF008F", "#FF0078", "#FF0060", "#FF0048", "#FF0030",
              "#FF0018")
    
  } else {
    
    
    if (colo=="fancy"){
      fancy.blue  <- c(c(255:0),rep(0,length(c(255:0))),rep(0,length(c(255:150))))
      fancy.green  <- c(c(0:255),c(255:0),rep(0,length(c(255:150))))
      fancy.red  <- c(c(0:255),rep(255,length(c(255:0))),c(255:150))
      colo <- rgb(b=fancy.blue/255,g=fancy.green/255,r=fancy.red/255)
      
    }
  }
  colorseq <- seq(0,1,length=length(colo))
  
  
  for (j in 1:dim(cl[[1]])[[1]]){
    if (single) j <- single
    tmp <- exprs(eset)[clusterindex==j, , drop=FALSE]# thanks Ian for the fix
    tmpmem <- memship[clusterindex==j,j]
    if (((j-1)%% (mfrow[1] * mfrow[2]))==0 | single){
      if (x11) X11(width=Xwidth,height=Xheight)
      #par(mfrow=mfrow)
      if (sum(clusterindex==j)==0) {
        ymin <- -1; ymax <- +1;
        
      } else {
        ymin <- min(tmp);ymax <- max(tmp);
        
      }
      
      
      if (sum(ylim.set == c(0,0)) ==2){
        ylim <- c(ymin,ymax)
      } else {
        ylim <- ylim.set
      }
      
      
      if (!is.na(sum(mfrow))){
        par(mfrow=mfrow,bg = bg,col.axis= col.axis,col.lab=col.lab,col.main=col.main,
            col.sub=col.sub,col=col)
      } else {
        par(bg = bg,col.axis= col.axis,col.lab=col.lab,col.main=col.main,
            col.sub=col.sub,col=col)
      }
      
      xlim.tmp <- c(1,dim(exprs(eset))[[2]])
      if (!(missing(time.points))) xlim.tmp <-  c(min(time.points), max(time.points))
      #print(c(j,ymin,ymax,ylim[1],ylim[2]));
      plot.default(x=NA,xlim=xlim.tmp, ylim= ylim,
                   xlab=xlab,ylab=ylab,main=paste("Cluster",j),axes=FALSE,...)
      
      if (missing(time.labels) && missing(time.points)){
        axis(1, 1:dim(exprs(eset))[[2]],c(1:dim(exprs(eset))[[2]]),col=ax.col,...)
        axis(2,col=ax.col,...)
      }
      if (missing(time.labels) && !(missing(time.points))) {
        axis(1, time.points,1:length(time.points),time.points,col=ax.col,...)
        axis(2,col=ax.col,...)
      }
      if (missing(time.points) & !(missing(time.labels))){
        axis(1, 1:dim(exprs(eset))[[2]],time.labels,col=ax.col,...)
        axis(2,col=ax.col,...)
      }
      if (!(missing(time.points)) & !(missing(time.labels))){  
        axis(1, time.points,time.labels,col=ax.col,...)
        axis(2,col=ax.col,...)
      }
      
      
    } else {
      
      if (sum(clusterindex==j)==0) {
        ymin <- -1; ymax <- +1;
        
      } else {
        ymin <- min(tmp);ymax <- max(tmp);
        
      }
      
      
      if (sum(ylim.set == c(0,0)) ==2){
        ylim <- c(ymin,ymax)
      } else {
        ylim <- ylim.set
      }
      
      xlim.tmp <- c(1,dim(exprs(eset))[[2]])
      if (!(missing(time.points))) xlim.tmp <-  c(min(time.points), max(time.points))
      #print(c(j,ymin,ymax,ylim[1],ylim[2]));  
      plot.default(x=NA,xlim=xlim.tmp, ylim= ylim,
                   xlab=xlab,ylab=ylab,main=paste("Cluster",j),axes=FALSE,...)
      
      
      if (missing(time.labels) && missing(time.points)){
        axis(1, 1:dim(exprs(eset))[[2]],c(1:dim(exprs(eset))[[2]]),col=ax.col,...)
        axis(2,col=ax.col,...)
      }
      if (missing(time.labels) && !(missing(time.points))) {
        axis(1, time.points,1:length(time.points),time.points,col=ax.col,...)
        axis(2,col=ax.col,...)
      }
      if (missing(time.points) & !(missing(time.labels))){
        axis(1, 1:dim(exprs(eset))[[2]],time.labels,col=ax.col,...)
        axis(2,col=ax.col,...)
      }
      if (!(missing(time.points)) & !(missing(time.labels))){  
        axis(1, time.points,time.labels,col=ax.col,...)
        axis(2,col=ax.col,...)
      }
      
      
    }
    
    if (length(tmpmem)>0){
      for (jj in 1:(length(colorseq)-1)){
        tmpcol <- (tmpmem >= colorseq[jj] & tmpmem <= colorseq[jj+1])
        
        if (sum(tmpcol)> 0) {
          tmpind <- which(tmpcol)
          for (k in 1:length(tmpind)){
            if (missing(time.points)) {
              lines(tmp[tmpind[k],],col=colo[jj],lwd=2)
            } else
              lines(time.points, tmp[tmpind[k],],col=colo[jj],lwd=2)
          }
        }
      }
    }
    #mat <- matrix(1:2,ncol=2,nrow=1,byrow=TRUE)
    #l   <- layout(mat,width=c(5,1))
    
    
    if (centre){
      
      lines(cl[[1]][j,],col=centre.col,lwd=centre.lwd)
      
    }
    
    
    #mfuzzColorBar()
    if (single) return();
    
  }
  
  
}



