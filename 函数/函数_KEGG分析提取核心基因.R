
prepare<- function() {
  project=""
  DEG=""
  n="" 
  hu_man=""
  human_db=""
  
  
}  


#默认为human，其他需要谨慎
KEGG_circ<- function(project,DEG,n, hu_man, human_db) {
  

    
    if (missing(hu_man)) 
      hu_man <- "human"
    
    if (missing(human_db)) 
      human_db <- "org.Hs.eg.db"
  
  
  if (!requireNamespace("ReactomePA", quietly = TRUE))
    BiocManager::install("ReactomePA",dependencies = TRUE)
  if (!requireNamespace(human_db, quietly = TRUE)) 
    BiocManager::install(human_db,dependencies = TRUE)
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) 
    BiocManager::install("clusterProfiler",dependencies = TRUE)
  if (!requireNamespace("enrichplot", quietly = TRUE)) 
    BiocManager::install("enrichplot",dependencies = TRUE)
  if (!requireNamespace("biomaRt", quietly = TRUE)) 
    BiocManager::install("biomaRt",dependencies = TRUE)
  if (!requireNamespace("tidyverse", quietly = TRUE)) 
    BiocManager::install("tidyverse",dependencies = TRUE)
  #data.table 更适用于处理大型数据集和需要高效率运算的情况，由于矛盾这两个只能liberary一个
  #而 dplyr 则更适合一般性的数据操作和数据清洗，尤其是对于需要清晰易读的代码和初学者而言。选择使用哪个取决于您的具体需求和个人偏好。
  library(tidyverse) 
  library(ReactomePA)
  library(enrichplot)
  library(biomaRt)
  library(org.Hs.eg.db) 
  library(clusterProfiler)
  library(dplyr)
  library(RColorBrewer) 
  
    library(GOplot)
    library(tidyr)
  # 1.数据处理（确定物种，更改OrgDb=org.Hs.eg.db。更改或检查DEG是否符合要求colnames(DEG)[1] <- "SYMBOL"，colnames(DEG)[2] <- "logFC"）
  # 重命名第一列
  colnames(DEG)[1] <- "SYMBOL"
  colnames(DEG)[2] <- "logFC"
  # 从第一列提取基因名称
  genename=as.character(DEG[, 1])
  gene_map=bitr(genename,fromType="SYMBOL",toType="ENTREZID",OrgDb=human_db)
  # 通过符号将gene_map与DEG合并
  gene_map <- inner_join(gene_map, DEG, by = "SYMBOL")
  # 移除含有NA值的行
  gene_map <- na.omit(gene_map)
  gene_map$SYMBOL <- toupper(gene_map$SYMBOL)

  # 获取排序后的索引
  sorted_index <- order(gene_map[, 3], decreasing = TRUE)
  # 根据排序后的索引重新排列数据框
  

  gene_map <- gene_map[sorted_index, ]
  # 提取logFC值，并将行名设为Entrez ID
  geneList <- gene_map[, 3]
  names(geneList) <- as.character(gene_map[, 2])
  
  
  
  
  # 2. 计算富集结果
  #计算exponent权重调节的是峰值的偏移一定程度会影响平滑程度，nPerm重复次数，大一点更平滑，但耗时

  # 2. 计算富集结果
  #2. 2 进行KEGG通路的GSEA
  KEGG_gseresult <- gseKEGG(geneList,organism = hu_man,exponent = 1.5, nPerm = 1000,maxGSSize = 1000, pvalueCutoff =1)
  KEGG_gseresult <- setReadable(KEGG_gseresult ,OrgDb=human_db,keyType ='ENTREZID')
  KEGG_gseresult <<-KEGG_gseresult
  Sys.sleep(5)
  
  KEGG_gseresult@result<-KEGG_gseresult@result[order(abs(KEGG_gseresult@result$enrichmentScore),decreasing = TRUE), ]
  
  Go1<-KEGG_gseresult@result[1:n,c(1,1,2,7,11)]
  #Go1<-Go_Reactomeresult@result[1:100,c(1,1,2,8,11)]
  colnames(Go1)<-c("category","ID","term","adj_pval","genes")
  # 假设您的数据框为 data，genes 列为需要修改的列
  Go1$genes <- gsub("/", ",", Go1$genes)
  #genes	A data frame with columns for 'ID', 'logFC'
  Go2<-gene_map[c(1,3)]#注意的是geneID是否是一样的格式，
  colnames(Go2)<-c("ID", "logFC")
  #合并获取最重要的数据框
  circ2 <-circle_dat(Go1,Go2)

 
  return(circ2)
}

#调用函数
#source("函数\\函数_KEGG分析提取核心基因.R")#默认为human，其他需要谨慎
#circ2=KEGG_circ(project,DEG,n=10, hu_man, human_db)