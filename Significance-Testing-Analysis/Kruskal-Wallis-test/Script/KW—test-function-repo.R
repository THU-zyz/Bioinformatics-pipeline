###### KW test function repo #####
library(readxl)
library(tidyverse)
library(dunn.test)

# Dataloader function 
DataLoader<-function(files,sheetnames){
  # Data preprocessing
  Data<-read_excel(files,sheet=sheetnames)
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

KW_test<-function(files,sheetnames,output){
  DataNclassnumber<-DataLoader(files,sheetnames)
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
  KW.Dunn<-KW.Dunn%>%select(`p-value: all`,`FDR: all`,everything())
  
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
  write.csv(KW.Dunn,output)
  return(KW.Dunn)
}