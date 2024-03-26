clean_data <- function(data) {
  # 查看数据结构
  str(data)
  
  # 检查是否有缺失值
  if (any(is.na(data))) {
    # 删除含有缺失值的行
    data <- data[complete.cases(data), ]
    
    # 删除含有缺失值的列
    data <- data[, colSums(is.na(data)) == 0]
    
    cat("含有缺失值的行/列已经被删除。\n")
  } else {
    cat("数据中没有缺失值。\n")
  }
  
  # 检查是否有重复的基因
  duplicated_genes <- rownames(data)[duplicated(rownames(data))]
  if (length(duplicated_genes) > 0) {
    print(paste("Duplicate genes found:", duplicated_genes))
    
    # 根据需要处理重复的基因
    for (gene_name in duplicated_genes) {
      # 找到重复的基因所在行的索引
      duplicated_rows <- which(rownames(data) == gene_name)
      
      # 计算重复行的平均值
      avg_values <- rowMeans(data[duplicated_rows, ])
      
      # 将平均值更新到第一个重复行
      data[duplicated_rows[1], ] <- avg_values
      
      # 删除其他重复行
      data <- data[-c(duplicated_rows[-1]), ]
    }
    cat("重复的基因已经被处理，取平均值。\n")
  } else {
    cat("数据中没有重复的基因。\n")
  }
  
  # log2转化
  if (max(data) > 50) {
    data <- log2(data + 1)
    cat("数据已经进行了log2转换。\n")
  }
  
  # 返回清洗后的数据
  return(data)
}
