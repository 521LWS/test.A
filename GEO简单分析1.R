# 如果没有安装BiocManager包，则安装它
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# 使用BiocManager安装所需的Bioconductor包
BiocManager::install("WGCNA", dependencies = TRUE)
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("affy")
BiocManager::install("tinyarray", dependencies = TRUE)
BiocManager::install("clusterProfiler")
BiocManager::install("ComplexHeatmap")
# 使用install.packages()函数安装其他需要的包
BiocManager::install("org.Rn.eg.db")
install.packages("AnnoProbe", dependencies = TRUE)
# 加载所需的R包
library("WGCNA")
library(tibble)
library("tinyarray")
library(stringr)
library(data.table)
library(AnnoProbe)
library(GEOquery)




# 下载并读取GEO数据集GSE35958
gse = "GSE35958"
geo = geo_download(gse,destdir=tempdir())
# 对表达数据进行log2转换（根据情况选择）
geo$exp <- log2(geo$exp + 1)

# 在Bioconductor中查找并加载GPL注释信息
find_anno(geo$gpl)
ids_570  <- AnnoProbe::idmap(geo$gpl,destdir = tempdir())




# 根据条件组（control和treatment）创建分组向量
group <- ifelse(str_detect(geo$pd$title, "ost"), "ost", "con")
# group 变量将被转换为一个因子，并且可以在后续的数据分析中使用
group <- factor(group)

# 提取符号（symbol）
# 当您需要处理和分析结构化数据，并且需要进行各种数据操作和统计分析时，
# 数据框是一个非常有用的数据结构。
ids <- as.data.frame(ids_570$SYMBOL)
# 当您需要存储和处理单一类型的一维数据集合，
# 并且进行数学运算、索引、切片、循环迭代等操作时，向量是一个非常有用的数据结构。
ids_v2 <- ids_570$SYMBOL

# 使用get_deg_all()函数找到所有的差异表达基因（DEGs）
DEGS  <-  get_deg_all(
  geo$exp,
  group,
  ids_570,
  symmetry = TRUE,
  my_genes = NULL,
  show_rownames = TRUE,
  cluster_cols = TRUE,
  color_volcano = c("#2874C5", "grey", "#f87669"),
  logFC_cutoff = 2,
  pvalue_cutoff = 0.05,
  adjust = FALSE,
  entriz = TRUE,
  n_cutoff = 2,
  annotation_legend = FALSE,
  lab = NA,
  species = "human"
)

# 生成差异表达基因的图像
DEGS$plots





