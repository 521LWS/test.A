# 1.准备工作
#通过禁止自动将字符向量转换为因子，你可以更好地控制数据处理过程，确保数据处理和建模的准确性和一致性。
options(stringsAsFactors = F)
# 如果没有安装BiocManager包，则安装它
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# 使用BiocManager安装所需的Bioconductor包
if (!requireNamespace("plyr", quietly = TRUE)) 
  BiocManager::install("plyr",dependencies = TRUE)
if (!requireNamespace("WGCNA", quietly = TRUE)) 
  BiocManager::install("WGCNA",dependencies = TRUE)
if (!requireNamespace("GEOquery", quietly = TRUE)) 
  BiocManager::install("GEOquery",dependencies = TRUE)
if (!requireNamespace("limma", quietly = TRUE)) 
  BiocManager::install("limma",dependencies = TRUE)
if (!requireNamespace("tinyarray", quietly = TRUE)) 
  BiocManager::install("tinyarray",dependencies = TRUE)
if (!requireNamespace("limma", quietly = TRUE)) 
  BiocManager::install("limma",dependencies = TRUE)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) 
  BiocManager::install("clusterProfiler",dependencies = TRUE)
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) 
  BiocManager::install("ComplexHeatmap",dependencies = TRUE)
if (!requireNamespace("org.Rn.eg.db", quietly = TRUE)) 
  BiocManager::install("org.Rn.eg.db",dependencies = TRUE)






# 使用install.packages()函数安装其他需要的包
if (!requireNamespace("AnnoProbe", quietly = TRUE)) 
     install.packages("AnnoProbe", dependencies = TRUE)
# 加载所需的R包

library(GEOquery)
library(tibble)
library("tinyarray")
library(AnnoProbe)

library(stringr)
library(data.table)

library("WGCNA")

# 2.下载，需要改gse
gse = "GSE35958"
geo = geo_download(gse,destdir=tempdir())
saveRDS(geo,file = paste0(gse,".rds"))
dim(geo$exp)
max(geo$exp)
# 根据max(geo$exp)的大小对表达数据进行log2转换（根据情况选择）
if(max(geo$exp)>50){#判断是否已经经过1og2转化
  geo$exp <- log2(geo$exp + 1)}#表达矩阵行1og2+1归一化处理
# 在Bioconductor中查找并加载GPL注释信息
checkGPL(geo$gpl)
printGPLInfo(geo$gpl)
find_anno(geo$gpl)
ids<- AnnoProbe::idmap(geo$gpl,destdir = tempdir())



# 3.分组，需要改分组关键词
group_list <- ifelse(str_detect(geo$pd$title, "ost"), "osteopo", "old")
# group 变量将被转换为一个因子，并且可以在后续的数据分析中使用
group_list <- factor(group_list)


# 4.作图，需要改参数
# 使用get_deg_all()函数找到所有的差异表达基因（DEGs）
DEGS  <-  get_deg_all(
  geo$exp,
  group_list,
  ids,
  
  symmetry = TRUE,
  my_genes = NULL,
  show_rownames = TRUE,
  cluster_cols = TRUE,
  color_volcano = c("#2874C5", "grey", "#f87669"),
  logFC_cutoff = 1,
  pvalue_cutoff = 0.05,
  adjust = FALSE,
  entriz = TRUE,
  n_cutoff = 2,
  annotation_legend =TRUE,
  lab = NA,
  species = "human"
)
# 生成差异表达基因的图像

filename <- gse
index <- 1
while (file.exists(paste0(filename, "-", index, ".pdf"))) {
  index <- index + 1
}
pdf(paste0(filename, "-", index, ".pdf"), onefile =TRUE,width =15, height =7.50)

DEGS$plots
dev.off()










