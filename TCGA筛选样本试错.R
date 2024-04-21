PROJECT <-'TCGA-LIHC'
counts<-`TCGA-LIHC_counts` #含有type，gene，未经处理
pd<-`TCGA-LIHC_pd`
h=150
keygene=c("PCK1","ODC1")
p1=as.data.frame(c("0","1"))
p2=as.data.frame(c("0","1"))
hu_man <- "human"
human_db <- "org.Hs.eg.db"



#1. 准备工作
options(scipen = 999)
# 获取当前代码文件的路径
current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# 将工作路径设置为当前代码文件所在的路径
setwd(current_dir)

library(data.table)
library(DESeq2)
library(tidyverse)
library(plyr)
library(tinyarray)

# 提取LUAD数据中的mRNA表达信息
counts_m <- counts[counts$gene_type == "protein_coding", ]
counts_m<- ddply(counts_m, .(gene_name), function(x) colMeans(x[, 3:ncol(x)]))
# 将第一列转换为行名
rownames(counts_m)<- counts_m$gene_name
# 删除原始的ID列
counts_m <- counts_m[, -1]
# 四舍五入表达矩阵中的值
counts_m<- round(counts_m, 0)



# # # # # # 分组信息
# 将mRNA表达数据的行名称转换为数据框形式
ids_m <- as.data.frame(rownames(counts_m))[c(1,1)]
# 根据mRNA表达数据创建TCGA组
group_list <- make_tcga_group(counts_m)
# 确定 "tumor" 和 "normal" 样本的列索引
tumor_index <- which(group_list  == "tumor")
normal_index <- which(group_list  == "normal")
# 重新排列数据框，确保 "tumor" 和 "normal" 样本连续出现
counts_m <- counts_m[, c(tumor_index, normal_index)]
# 根据mRNA表达数据创建TCGA组
group_list <- make_tcga_group(counts_m)

colData <-as.data.frame(colnames(counts_m)) 
colData$Type<-group_list 
# 修改第二列名为 'Type'
colnames(colData)[c(1,2)] <- c("sample","Type")
# 确认分组信息中的水平


# 2. 标化处理(关键在于更改样品筛选中的h值)
# 创建 DESeq2 数据集对象
dds <- DESeqDataSetFromMatrix(countData = counts_m, colData = colData, design = ~ Type)
# 过滤低表达的基因

dds <- dds[rowSums(counts(dds)) > 1,]
# VST 标准化处理

vsd <- vst(dds, blind = FALSE)
expr_matrix  <-assay(vsd)

# 使用 plotPCA 进行主成分分析
plotPCA(vsd, "Type")












# 定义 h 的范围
h_values <- seq(116, 124, by = 1)

# 遍历每个 h 值并运行程序
for (i in 1:length(h_values)) {
  # 在此处运行您的程序，使用当前的 h 值
  h=h_values[i]

# 2.1 样品筛选(关键在于更改样品筛选中的h值)

source("函数\\函数_通过聚类进行样品筛选.R")
process_sample_data<-process_sample_data(expr_matrix,h)
processed_counts_m<-process_sample_data


#3. 重新分组计算
# # # # # # 分组信息
# 将mRNA表达数据的行名称转换为数据框形式
ids_m <- as.data.frame(rownames(processed_counts_m))[c(1,1)]
# 根据mRNA表达数据创建TCGA组
group_list <- make_tcga_group(processed_counts_m)
# 确定 "tumor" 和 "normal" 样本的列索引
tumor_index <- which(group_list  == "tumor")
normal_index <- which(group_list  == "normal")
# 重新排列数据框，确保 "tumor" 和 "normal" 样本连续出现
processed_counts_m <- processed_counts_m[, c(tumor_index, normal_index)]
processed_counts_m <-as.data.frame(processed_counts_m)
# 根据mRNA表达数据创建TCGA组
group_list <- make_tcga_group(processed_counts_m)
print(group_list)
colData <-as.data.frame(colnames(processed_counts_m)) 
colData$Type<-group_list 
# 修改第二列名为 'Type'
colnames(colData)[c(1,2)] <- c("sample","Type")
# 确认分组信息中的水平


row_indices <- match(colnames(processed_counts_m),colnames(counts_m))
# 根据表B的行索引对表A的行进行排序
counts_m1<-counts_m[, row_indices]
dds <- DESeqDataSetFromMatrix(countData = counts_m1, colData = colData, design = ~ Type)
# 过滤低表达的基因
dds <- dds[rowSums(counts(dds)) > 1,]
# VST 标准化处理
vsd <- vst(dds, blind = FALSE)
# 使用 plotPCA 进行主成分分析
plotPCA(vsd, "Type")


# 进行差异分析(耗时较长)
dds2 <- DESeq(dds)
# 获取差异分析结果
res <- results(dds2)
# 设置对照组
contrast <- c("Type", "tumor", "normal")
# 获取差异分析结果
DEG <- as.data.frame(results(dds2, contrast = contrast, alpha = 0.05))
# 将结果转换为 tibble 格式(行名入表命名为gene)
DEG <- tibble::rownames_to_column(DEG, var = "Gene")
DEG1<-DEG[c(1,3,7)]
colnames(DEG1)<-c("Gene","logFC","pvalue")
# 导出差异分析结果
df<- subset(DEG1,abs(logFC) >= 0 )
# 例如：h = h

p1[1, i] <- h_values[i]
for (n in 1:length(keygene)) {
  cat("Looking for", keygene[n], "\n")
  matches <- df[[1]] == keygene[n]
  if (any(matches)) {
    p1[n+1, i] <- df[matches, 2]
    cat("Found:", df[matches, 2], "\n")
  } else {
    cat("No match found for", keygene[n], "\n")
  }
  
  
  p2[1, i] <- h_values[i]
  for (n in 1:length(keygene)) {
    cat("Looking for", keygene[n], "\n")
    matches <- df[[1]] == keygene[n]
    if (any(matches)) {
      p2[n+1, i] <- df[matches, 3]
      cat("Found:", df[matches, 3], "\n")
    } else {
      cat("No match found for", keygene[n], "\n")
    }
  
  
}



saveRDS(p1, file = paste0("h-","p1.rds"))

saveRDS(df, file = paste0(h,"df.rds"))

}