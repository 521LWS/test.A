
prepare<- function() {
  geo_data=""
  h=""
  
  
}  

process_sample_data <- function(geo_data, h) {
  library("WGCNA")

  library(AnnoProbe)
  # 根据MAD方法筛选变化大的基因并转置
  processed_geo <- t(geo_data[order(apply(geo_data, 1, mad), decreasing = TRUE)[1:10000], ])
  
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
  sampleTree<<-  hclust(dist(processed_geo), method = "average")
  sampleTree<-  hclust(dist(processed_geo), method = "average")
  par(pin = c(12, 9), cex = 0.6, mar = c(0, 4, 2, 0))
  plot(sampleTree, main = "Sample clustering after removing outliers", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
  
  # 按照样本聚类结果对基因表达数据进行排序
  row_indices <- match(rownames(processed_geo), colnames(geo_data))
  processed_geo <- geo_data[, row_indices]
  return(processed_geo)

}
#调用函数(需要根据输出或报错调整h的大小)
#source("函数\\函数_通过聚类进行样品筛选.R")
#process_sample_data =process_sample_data(geo_data=counts_m,h=240)

