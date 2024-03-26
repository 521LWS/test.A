
#下载下载器TCGAbiolinks，当更新出现问题，如没有挪动的权限，直接点右侧手动更新即可
package_name  <- "BiocManager"
if (!requireNamespace(package_name , quietly = TRUE)) 
  install.packages(package_name )
package_name  <- "TCGAbiolinks"
if (!requireNamespace(package_name , quietly = TRUE)) 
  BiocManager::install(package_name )

library(TCGAbiolinks)
library(SummarizedExperiment)







getGDCprojects()$project_id 
CANCER <- GDCquery(project = "TCGA-COAD", 
                 data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "STAR - Counts")
#设置 method = "api" 是为了确保您从 GDC 获取数据时使用的是 API 方法，而不是其他方法。这对于确保数据下载过程的可靠性和数据的一致性非常重要。
GDCdownload(query = CANCER, method = "api")
#合并并创建一个对象储存数据#可以简单理解，我们前面的下载只不过是下载了一个个样本，然后通过上面这个函数把所有样本进行整合，变成了一个整体的对象
expr <- GDCprepare(query=CANCER)
#.Counts = "unstranded" tpm = "tpm_unstrand" fpkm = " fpkm_unstrand"
counts <- as.data.frame(assay(expr))#默认提取counts数据
TPM <- as.data.frame(assay(expr,i = "tpm_unstrand"))#提取TPM数据
data=as.data.frame(rowRanges(expr))#获取其它信息数据#这里面就包括注释以及编码、非编码等等信息（包括标准化所需的基因长度等）
expr_count = cbind(gene_type=data$gene_type,gene_name=data$gene_name,counts)
counts<-TPM

#数据清洗：data_cleaning.R文件
source("data_cleaning.R")
cleaned_data <- clean_data(counts)

# 加载所需的包
library(data.table)
library(tinyarray)

# 从文件中读取LUAD数据并转换为数据框
counts <- as.data.frame(counts)
# 提取LUAD数据中的lncRNA表达信息
counts_lnc <- trans_exp(counts, lncrna_only = TRUE)
# 提取LUAD数据中的mRNA表达信息
counts_m <- trans_exp(counts, mrna_only = TRUE)

# 根据mRNA表达数据创建TCGA组
group <- make_tcga_group(counts_m)
# 将mRNA表达数据的行名称转换为数据框形式
ids_m <- as.data.frame(rownames(counts_m))
ids_m <-ids_m [c(1,1)]
# 获取mRNA表达数据的行名称
ids_mSv2 <- rownames(counts_m)

degs_m <- get_deg_all(
  exp = counts_m,  # 使用LUAD_m作为表达数据
  group_list = group,  # 使用之前创建的group作为分组信息
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
  n_cutoff = 3,
  annotation_legend = FALSE,
  lab = NA,
  species = "human"
)

degs_m$plots


  