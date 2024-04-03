# 1.准备工作
#下载下载器TCGAbiolinks，当更新出现问题，如没有挪动的权限，清除环境数据后，直接点右侧手动更新即可
package_name  <- "BiocManager"
if (!requireNamespace(package_name , quietly = TRUE)) 
  install.packages(package_name )
package_name  <- "TCGAbiolinks"
if (!requireNamespace(package_name , quietly = TRUE)) 
  BiocManager::install(package_name )
if (!requireNamespace("data.table", quietly = TRUE)) 
  BiocManager::install("data.table",dependencies = TRUE)
if (!requireNamespace("tinyarray", quietly = TRUE)) 
  BiocManager::install("tinyarray",dependencies = TRUE)



library(TCGAbiolinks)
library(SummarizedExperiment)
projects<-getGDCprojects()$project_id 
print(projects)
projects <- projects[grepl('^TCGA', projects , perl=TRUE)]
print(projects)




# 2.下载临床数据和表达数据（只需要更改project）
#获取临床数据，新版用TCGA_pd <- GDCquery_clinic(project = "TCGA-COAD")，存在无效行
#参考??GDCprepare_clinic
project <- "TCGA-LUAD"
query <- GDCquery(project = project,
                  data.category = "Clinical", 
                  data.type = "Clinical Supplement",
                  data.format ="bcr xml")
GDCdownload(query)
#患者信息 (clinical.info = "patient")：如患者ID、年龄、性别、生存时间、肿瘤分期等。这些信息通常用于了解患者群体的特征，以及评估与生存率、治疗反应等相关的临床特征。
#药物信息 (clinical.info = "drug"):如药物名称、给药剂量、给药时间、药物反应等。这些信息对于了解患者的治疗历史以及评估不同治疗方案的效果非常重要。
#放射治疗信息 (clinical.info = "radiation")：如治疗剂量、治疗时间、治疗区域等。这些信息对于了解患者的治疗历史、评估放射治疗的效果以及预测患者的预后都非常重要。
#行政信息 (clinical.info = "admin")：如入院时间、出院时间、医院名称、病房号等。这些信息通常用于管理患者的治疗过程、跟踪患者的就诊情况以及管理临床试验的进行。
pd <- GDCprepare_clinic(query, clinical.info = "patient")


#获取表达数据
project = "TCGA-LUAD"
CANCER <- GDCquery(project = project, 
                 data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "STAR - Counts")
#设置 method = "api" 是为了确保您从 GDC 获取数据时使用的是 API 方法，而不是其他方法。这对于确保数据下载过程的可靠性和数据的一致性非常重要。
GDCdownload(query = CANCER, method = "api")
#合并并创建一个对象储存数据#可以简单理解，我们前面的下载只不过是下载了一个个样本，然后通过上面这个函数把所有样本进行整合，变成了一个整体的对象
expr <- GDCprepare(CANCER,save = T,save.filename = paste0(project,"_RNA.Rdata"))
#.Counts = "unstranded" tpm = "tpm_unstrand" fpkm = " fpkm_unstrand"
#TPM和FPKM没有考虑到高表达基因波动以及有效reads总数对基因表达量的影响。
expr=data
counts <- as.data.frame(assay(expr))#默认提取counts数据
TPM <- as.data.frame(assay(expr,i = "tpm_unstrand"))#提取TPM数据
data1=as.data.frame(rowRanges(expr))#获取其它信息数据#这里面就包括注释以及编码、非编码等等信息（包括标准化所需的基因长度等）
expr_count = cbind(gene_type=data1$gene_type,gene_name=data1$gene_name,counts)

# 3.转换和处理数据格式(需要选择是否清洗数据，以及确认清洗时基因列的列名是否正确)
# 加载所需的包
library(data.table)
library(tinyarray)

#清洗处理数据并
counts<-expr_count
#source("函数\\函数_data_cleaning.R")
#counts<- clean_data(counts, 'gene_name')



# 提取LUAD数据中的lncRNA表达
counts_lnc <- expr_count [expr_count$gene_type == "lncRNA", ]
counts_lnc <- counts_lnc[, -1]
#数据清洗：data_cleaning.R文件
#source("函数\\函数_data_cleaning.R")
#counts_lnc_cleaning<- clean_data(counts_lnc, 'gene_name')
# 提取LUAD数据中的mRNA表达信息
counts_m <- expr_count [expr_count$gene_type == "protein_coding", ]
counts_m <- counts_m[, -1]




# 4.分组（在此处输入已处理过的表达矩阵）
#counts_m <-`TCGA-COAD_counts_m`
#数据清洗：data_cleaning.R文件
source("函数\\函数_data_cleaning.R")
counts_m_cleaning<- clean_data(counts_m, 'gene_name')
# 将mRNA表达数据的行名称转换为数据框形式
ids_m <- as.data.frame(rownames(counts_m_cleaning))[c(1,1)]

# 根据mRNA表达数据创建TCGA组
group_list <- make_tcga_group(counts_m_cleaning)
# 确定 "tumor" 和 "normal" 样本的列索引
tumor_index <- which(group_list  == "tumor")
normal_index <- which(group_list  == "normal")
# 重新排列数据框，确保 "tumor" 和 "normal" 样本连续出现
counts_m_cleaning <- counts_m_cleaning[, c(tumor_index, normal_index)]
# 根据mRNA表达数据创建TCGA组
group_list <- make_tcga_group(counts_m_cleaning)


# 5.作图
DEGS <- get_deg_all(
  exp = counts_m_cleaning,  # 使用LUAD_m作为表达数据
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

# 6. 储存
# 创建一个名为 "geo" 的文件夹
dir.create(project)
saveRDS(counts,file = paste0(project,"/",project,"_counts.rds"))
saveRDS(counts_m,file = paste0(project,"/",project,"_counts_m.rds"))
saveRDS(counts_lnc,file = paste0(project,"/",project,"_counts_lnc.rds")) 
saveRDS(counts_m_cleaning,file = paste0(project,"/",project,"_counts_m_cleaning.rds")) 
saveRDS(pd,file = paste0(project,"/",project,"_pd.rds")) 

filename <- project
index <- 1
while (file.exists(paste0(filename, "-", index, ".pdf"))) {
  index <- index + 1
}
pdf(paste0(filename,"/",filename, "-", index, ".pdf"), onefile =TRUE,width =15, height =7.50)

DEGS$plots
dev.off()


