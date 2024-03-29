
process_sample_data <- function(geo_data, h) {
  library("WGCNA")

  library(AnnoProbe)
  source("data_cleaning.R")
  
  # 清理基因表达数据
  processed_geo <- clean_data(geo_data)
  
  # 根据MAD方法筛选变化大的基因并转置
  processed_geo <- t(processed_geo[order(apply(processed_geo, 1, mad), decreasing = TRUE)[1:10000], ])
  
  # 对样本进行聚类以检测异常值
  sampleTree <- hclust(dist(processed_geo), method = "average")
  par(pin = c(12, 9), cex = 0.6, mar = c(0, 4, 2, 0))
  plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  abline(h = h, col = "red")
  
  # 根据聚类结果筛选样本
  keepSamples <- (cutreeStatic(sampleTree, cutHeight = h, minSize = 4) == 1)
  print(keepSamples)
  processed_geo <- processed_geo[keepSamples, ]
  
  # 重新计算样本聚类树
  sampleTree <- hclust(dist(processed_geo), method = "average")
  par(pin = c(12, 9), cex = 0.6, mar = c(0, 4, 2, 0))
  plot(sampleTree, main = "Sample clustering after removing outliers", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  
  # 按照样本聚类结果对基因表达数据进行排序
  row_indices <- match(rownames(processed_geo), colnames(geo_data))
  processed_geo <- geo_data[, row_indices]
  processed_geo <<- processed_geo
  return(NULL)

}



