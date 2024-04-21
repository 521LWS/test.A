

prepare<- function() {
  data=""
  gene="GENE"
  
  
}  



clean_data <- function(data,gene) {
  # 查看数据结构
  str(data)
  
  # 检查是否有缺失值
  if (any(is.na(data))) {
    # 删除含有缺失值的行
   data <- na.omit(data)  
   data <- data[complete.cases(data), ]
   #删除整行都是缺失值的行
   data <- data[rowSums(is.na(data)) != ncol(data), ]

  
    # 删除含有缺失值的列
    data <- data[, colSums(is.na(data)) == 0]
    
    cat("含有缺失值的行/列已经被删除。\n")
  } else {
    cat("数据中没有缺失值。\n")
  }
  
  
  # 检查是否存在 "gene" 列
  if (gene %in% colnames(data)) {
    # 检查 "gene" 列中是否有重复的值
    duplicated_values <- duplicated(data[, gene])
    gene_index <- grep(gene, colnames(data)) 
    colnames(data)[gene_index] <- "gene"
    if (any(duplicated_values)) {
      print(paste("Duplicate values found in the gene column:", data[duplicated_values, gene]))
      
      # 根据需要处理重复的值
      library(plyr)
      # 按基因符号计算平均表达量
      data <- ddply(data, .(gene), function(x) colMeans(x[, (gene_index+1):ncol(x)]))
      # 将第一列转换为行名
      rownames(data) <- data$gene
      data <- data[, -1]
      cat("重复的值已经被处理。\n")
    } else {
      rownames(data) <- data$gene
      data <- data[, -gene_index:-1]
      cat("数据中没有重复的值。\n")
    }
  } else {
    cat("数据中不存在 gene 列。\n")
  }
  
  
      
  # 删除所有元素都小于等于0的行
  if (any(data <= 0)) {
    #rowSums(data <= 0) != ncol(data)：这部分代码与 ncol(data)（数据框的列数）进行比较，以确定哪些行中至少有一个元素小于等于 0。
    data <- data[rowSums(data <= 0) != ncol(data), ]
    cat("数据已经删除小于等于0的行。\n")
  }
  
  
  # log2转化

  if (max(data) > 300) {
    data <- log2(data + 1)
    cat("数据已经进行了log2转换。\n")
  }
  # 返回清洗后的数据
  return(data)
}
#调用函数
#source("函数\\函数_data_cleaning.R")#去除NA，去除gene列的重复行取平均，取log
#processed_geo <- clean_data(data=geo_exp, gene="symbol")