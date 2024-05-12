

prepare<- function() {
  geo_exp=""
  ids=""
  
  
}  

process_geo <- function(geo_exp,ids) {
  library("tinyarray")
  library(AnnoProbe)
  library(plyr)
 
  
  # 提取基因表达数据
  geo_exp <- as.data.frame(geo_exp)
  geo_exp <- tibble::rownames_to_column(geo_exp, var = "Gene")
  
  # 合并基因表达数据和探针注释
  geo_exp <- merge(ids, geo_exp, by.x = "probe_id", by.y = "Gene")
  
  # 按基因符号计算平均表达量
  geo_exp <- ddply(geo_exp, .(symbol), function(x) colMeans(x[, 3:ncol(x)]))
  
  # 将第一列转换为行名
  rownames(geo_exp) <- geo_exp$symbol
  
  # 删除原始的ID列
  geo_exp <- geo_exp[, -1]
  
  return(geo_exp)
}

# 调用函数
# source("函数\\函数_geo探针替换.R")#要求gene栏入行名,ids列名为probe_id,symbol
# processed_geo <- process_geo(geo_exp,ids)

