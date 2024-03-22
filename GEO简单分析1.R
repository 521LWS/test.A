
#通过禁止自动将字符向量转换为因子，你可以更好地控制数据处理过程，确保数据处理和建模的准确性和一致性。
options(stringsAsFactors = F)
# 如果没有安装BiocManager包，则安装它
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
# 使用BiocManager安装所需的Bioconductor包
BiocManager::install("plyr")
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
library(plyr)


getoption('timeout')
options(timeout=100000)
# 下载并读取GEO数据集GSE35958



gse = "GSE35958"
geo = geo_download(gse,destdir=tempdir())
dim(geo$exp)
max(geo$exp)
# 根据max(geo$exp)的大小对表达数据进行log2转换（根据情况选择）
if(max(geo$exp)>50){#判断是否已经经过1og2转化
  geo$exp <- log2(geo$exp + 1)#表达矩阵行1og2+1归一化处理
# 在Bioconductor中查找并加载GPL注释信息
checkGPL(geo$gpl)
printGPLInfo(geo$gpl)
find_anno(geo$gpl)
ids_570<- AnnoProbe::idmap(geo$gpl,destdir = tempdir())
ids <- as.data.frame(ids_570$SYMBOL)



group <- ifelse(str_detect(geo$pd$title, "ost"), "ost", "con")
# group 变量将被转换为一个因子，并且可以在后续的数据分析中使用
group <- factor(group)



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
  logFC_cutoff = 1,
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















# 数据处理
geo_exp <- geo$exp
geo_pd <- geo$pd
# 更改格式并命名
geo_exp <- as.data.frame(geo_exp)
geo_exp <- tibble::rownames_to_column(geo_exp, var = "Gene")
# 合并基因表达数据和探针注释
geo_exp <- merge(ids_570, geo_exp, by.x = "probe_id",by.y ="Gene")
# 按基因符号计算平均表达量
geo_exp <- ddply(geo_exp, .(symbol), function(x) colMeans(x[, 3:ncol(x)]))



# 分组：提取第一列并新建group_list.总体不如group <- ifelse(str_detect(geo$pd$title, "ost"), "ost", "con")# group 变量将被转换为一个因子，并且可以在后续的数据分析中使用group <- factor(group)
group_list=geo_pd [1]
colnames(group_list)="group"
# 三种文字编辑方式根据需要选择，其实不如上面的判断句来的好用
group_list$group <- stringr::str_remove(group_list$group, "hMSC")
group_list$group <- str_replace(group_list$group, ".*-(.*?)_.*", "\\1")
group_list$group <- str_extract(group_list$group, ".*-(.*?)_")
#group_list <- group_list[colnames(geo_exp),]：这种写法中，逗号 , 前面的部分指定了要选择的行，而逗号后面的部分为空，这意味着你选择了所有的列。这样做会保留 group_list 中与 geo_exp 中列名相匹配的行，并且保留所有的列。
#group_list <- group_list[colnames(geo_exp)]：这种写法中，逗号后面的部分被省略了，这意味着你只选择了列，而没有选择行。这样做会保留 group_list 中与 geo_exp 中列名相匹配的列，并且保留所有的行。
#总之，区别在于逗号后面是否带有空值，决定了你是选择了行还是列。其实这个函数意义不大，更严谨。
#这个操作可能会导致 group_list 中的行名和顺序与 geo_exp 中的列名和顺序相匹配从而达到正确调用sample的目的。
group_list <-group_list[colnames(geo_exp),]
table(group_list)
print(group_list)


# 加载所需的包
library(limma)
# 创建设计矩阵（更改分组group_list）
design <- model.matrix(~0 + factor(group_list))
colnames(design) <- levels(factor(group_list))




# 设置对比矩阵,control在后
formula <- "osteopo-old"
contrast.matrix <- makeContrasts(contrasts= formula, levels = design)
print(contrast.matrix)
# 拟合线性模型
fit <- lmFit(geo_exp, design)
# 进行对比
fit2 <- contrasts.fit(fit, contrast.matrix)
# 应用贝叶斯估计
fit2 <- eBayes(fit2)
# 获取差异表达基因（DEG）
DEG <- topTable(fit2, coef = 1, n = Inf)
DEG <- na.omit(DEG)
# 将结果转换为 tibble 格式(行名入表命名为gene)
DEG <- tibble::rownames_to_column(DEG, var = "Gene")
# 整理并重命名
df<-DEG[c(1,2,5)]
colnames(df)<-c("Gene","10g2FC","pvalue")




