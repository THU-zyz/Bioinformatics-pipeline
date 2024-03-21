library(rio)
library(clusterProfiler)

setwd("~/Desktop/博士生活/科研生活/Hu_lab/Github/Bioinformatics-pipeline/Pathway-analysis/KEGG-Pathway-Enrichment/Metabolomics/")
##### Multi-Group KEGGEnrich #####
MetabolicPathway<-import("Data/MetabolicPathway.xlsx")
path.cpd<-import("Data/path.cpd.xlsx")

MultiGroupKEGGEnrich<-function(Cluster,MetabolicPathway,path.cpd){
    
}
