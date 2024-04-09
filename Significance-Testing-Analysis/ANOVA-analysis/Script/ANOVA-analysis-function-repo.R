##### ANOVA analysis #####
library(tidyverse)
ANOVA_analysis<-function(file_path,output_name){
  Data<-read.xlsx(file_path,sheetIndex = 1)
  # 确保group列是因子类型
  Data$group<-as.factor(Data$group)
  # 初始化一个空的数据框来收集ANOVA结果
  anova_results <- data.frame()
  # 对每个代谢物进行ANOVA分析
  for (metabolite in names(Data)[3:ncol(Data)]) {
    
    result <- aov(Data[[metabolite]] ~ group, data = Data)
    anova_table <- summary(result)[[1]]
    
    # 将结果添加到anova_results数据框
    anova_results <- rbind(anova_results, 
                           cbind(Metabolite = metabolite, 
                                 anova_table[1, ], 
                                 Pr_Gt_F = anova_table[1, "Pr(>F)"]))
  }
  anova_results<-anova_results%>%dplyr::select(Metabolite,SumSq=`Sum Sq`,MeanSq=`Mean Sq`,Fvalue=`F value`,Pvalue=`Pr(>F)`)
  rownames(anova_results) <- NULL
  anova_results$p.adjust<-p.adjust(anova_results$Pvalue,method = "BH")
  anova_results<- anova_results[order(anova_results$p.adjust), ]
  write.xlsx(anova_results,file = paste0("Result/",output_name,"-ANOVA-result.xlsx"))
  print("ANOVA analysis finished successfully!")
}