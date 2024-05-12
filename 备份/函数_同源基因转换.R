if (!requireNamespace("biomaRt", quietly = TRUE))
  BiocManager::install("biomaRt",dependencies = TRUE)
library(biomaRt)
# 指定 Ensembl 镜像服务器
ensembl <- useEnsembl(biomart = "ensembl", mirror = "asia")
# 指定用于查询的人类和小鼠数据库。查了谷歌后发现报错是网页自身的问题，在构建数据集的时候需更换2021年版本的一个网页才能正常运行，估计是2022年版本的bug。
#解决办法是是用host参数指定2021年版本网页。
human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
rat_mart <- useMart("ensembl", dataset = "rnorvegicus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 


mouse_to_human <- function(mouse_genes) {

  # 获取小鼠基因到人类基因的映射信息
  mouse_to_human <- getLDS(
    attributes = c("mgi_symbol"),  # 小鼠基因属性
    filters = "mgi_symbol",  # 小鼠基因过滤器
    values =mouse_genes ,  # 小鼠基因列表
    mart = mouse_mart,  # 小鼠数据库
    attributesL = c("hgnc_symbol"),  # 人类基因属性
    martL =  human_mart,  # 人类数据库
    uniqueRows = TRUE
  )
  # 显示结果
  head(mouse_to_human)
  colnames(mouse_to_human)<-c("probe_id","symbol")
  return(mouse_to_human)
}

  
  
rat_to_human <- function(rat_genes) {

  # 获取大鼠基因到人类基因的映射信息
  rat_to_human <- getLDS(
    attributes = c("rgd_symbol"),  # 大鼠基因属性
    filters = "rgd_symbol",  # 大鼠基因过滤器
    values = rat_genes,  # 大鼠基因列表
    mart = rat_mart,  # 大鼠数据库
    attributesL = c("hgnc_symbol"),  # 人类基因属性
    martL = human_mart,  # 人类数据库
    uniqueRows = TRUE
  )
  # 显示结果
  head(rat_to_human)
  colnames(rat_to_human)<-c("probe_id","symbol")
  return(rat_to_human)
} 
  
  



#<- 
## 调用函数
# source("函数\\函数_同源基因转换.R")
# mouse_to_human <- mouse_to_human(mouse_genes=rownames(GSE126500$exp))

# rat_to_human <- rat_to_human(rat_genes=rownames(GSE126500$exp))

