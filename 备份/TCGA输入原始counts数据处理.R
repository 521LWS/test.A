#counts <- as.data.frame(assay(expr))#默认提取counts数据

打开TCGA-LUAD_pd，打开TCGA-LUAD_counts，先修改以下内容


PROJECT <-'TCGA-LUAD'
counts<-`TCGA-LUAD_counts`#含有type，gene，未经处理
pd<-`TCGA-LUAD_pd`
h=155




#1. 准备工作（需要修改counts_m）
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









# 2.1 样品筛选(关键在于更改样品筛选中的h值)

source("函数\\函数_通过聚类进行样品筛选.R")
process_sample_data =process_sample_data(expr_matrix,h)
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
# 根据mRNA表达数据创建TCGA组
group_list <- make_tcga_group(processed_counts_m)
print(group_list)
colData <-as.data.frame(colnames(processed_counts_m)) 
colData$Type<-group_list 
# 修改第二列名为 'Type'
colnames(colData)[c(1,2)] <- c("sample","Type")
# 确认分组信息中的水平



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
df<- subset(DEG1,abs(logFC) >= 0 & pvalue< 0.05)




#4. 处理临床信息（已删除无效值，但没有去除异常样品，还需后期处理）


source("函数\\函数_TCGA临床信息过滤.R")
process_TCGA_pd <- process_TCGA_clinical(counts_m,pd)
numeric_TCGA_pd <- numeric_TCGA_clinical(`process_TCGA_pd`)

# 获取表B的行索引
#row_indices <- match(colnames(processed_counts_m), rownames(process_TCGA_pd))
# 根据表B的行索引对表A的行进行排序
#processed2_pd <-process_TCGA_pd[row_indices, ]




#.6WGCNA分析
source("函数\\函数_WGCNA出图.R")
WGCNA(project=PROJECT, `processed_counts_m`, `numeric_TCGA_pd`,module="turquoise")



















#（此处可以暂停，module="turquoise"需要手动更改生成的select-gene需要手动筛选模块儿，可以利用该基因组进行df的筛选）
#df <- df[match(select_gene, df$Gene), ]
#df <- na.omit(df)


#.7.1聚类分析
source("函数\\函数_聚类分析.R")#默认为human，其他需要谨慎
go_KEGG(project=PROJECT,DEG=df)


#7.2 环形热图(利用df，geo_exp)(样品太多，这里只选了10个,注意超过10个由于默认按文字排序，会出现1,10,2的顺序)

processed_counts_m1=processed_counts_m[ , c(496:504)]

#调用函数
source("函数\\函数_环形热图_依据KEGG筛选的gene.R")#默认为human，其他需要谨慎
heatmap_KEGG(project=PROJECT,EXP=processed_counts_m1,DEG=df)

#7.3 免疫浸润分析（需要自己选基因）

select_gene=df[1:20,1]
expr=processed_counts_m
geneset = rio::import("Marker.xlsx")
Group<- tinyarray::make_tcga_group(expr)

#source("函数\\函数_免疫浸润分析.R")
#immune_cell(PROJECT,expr,geneset,Group,select_gene)





# 9.作图
DEGS <- get_deg_all(
  exp = processed_counts_m,  # 使用LUAD_m作为表达数据
  group_list = group_list,  # 使用之前创建的group作为分组信息
  ids = ids_m,  # 使用之前创建的ids_m作为ID信息
  logFC_cutoff = 1,
  pvalue_cutoff = 0.05,
  adjust = TRUE,
  entriz = FALSE,
  symmetry = FALSE,  # 将对称性设置为FALSE
  my_genes = NULL,
  show_rownames = FALSE,
  cluster_cols = FALSE,  # 将cluster_cols设置为FALSE
  color_volcano = c("#2874c5", "grey", "#f87669"),
  n_cutoff = 2,
  annotation_legend = FALSE,
  lab = NA,
  species = "human"
)

DEGS$plots



#.10.储存

# 创建一个名为 "geo" 的文件夹
dir.create(PROJECT)
# 保存 RDS 文件到新建的文件夹中
saveRDS(processed_counts_m, file = paste0(PROJECT,"/",PROJECT,"processed_counts_m.rds"))
saveRDS(DEG1, file = paste0(PROJECT,"/",PROJECT,"DEG1(mean).rds"))
saveRDS(df, file = paste0(PROJECT,"/",PROJECT,"df.rds"))

saveRDS(numeric_TCGA_pd, file = paste0(PROJECT,"/",PROJECT,"numeric_TCGA_pd.rds"))
saveRDS(process_TCGA_pd, file = paste0(PROJECT,"/",PROJECT,"process_TCGA_pd.rds"))

# 存图，需要改名字
#如果执行 pdf() 函数后输出的 PDF 文件是空的，可能是因为在调用 pdf() 函数之后没有进行任何绘图操作，或者绘图操作没有被保存到 PDF 文件中。
filename <- PROJECT
index <- 1
while (file.exists(paste0(filename, "-", index, ".pdf"))) {
  index <- index + 1
}
pdf(paste0(filename,"/",filename, "-", index, ".pdf"), onefile =TRUE,width =15, height =7.50)

DEGS$plots
dev.off()

