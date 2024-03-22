
#下载下载器TCGAbiolinks，当更新出现问题，如没有挪动的权限，直接点右侧手动更新即可
package_name  <- "BiocManager"
if (!requireNamespace(package_name , quietly = TRUE)) 
  install.packages(package_name )
package_name  <- "TCGAbiolinks"
if (!requireNamespace(package_name , quietly = TRUE)) 
  BiocManager::install(package_name )


library("TCGAbiolinks")
library(TCGAbiolinks)
library(SummarizedExperiment)

getGDCprojects()$project_id 
COAD <- GDCquery(project = "TCGA-COAD", 
                 data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification", 
                 workflow.type = "STAR - Counts")

GDCdownload(query = COAD, method = "api")
#合并并创建一个对象储存数据#可以简单理解，我们前面的下载只不过是下载了一个个样本，然后通过上面这个函数把所有样本进行整合，变成了一个整体的对象
expr <- GDCprepare(query=COAD)

#.Counts = "unstranded" tpm = "tpm_unstrand" fpkm = " fpkm_unstrand"
counts <- as.data.frame(assay(expr))#默认提取counts数据
TPM <- as.data.frame(assay(expr,i = "tpm_unstrand"))#提取TPM数据
data=as.data.frame(rowRanges(expr))#获取其它信息数据#这里面就包括注释以及编码、非编码等等信息（包括标准化所需的基因长度等）
expr_count = cbind(gene_type=data$gene_type,gene_name=data$gene_name,counts)





library(data.table)
library(tinyarray)
LUAD <-fread("TCGA-LUAD.htseq_counts.tsv.gz")
LUAD <-as.data.frame(LUAD)
LUAD<-LUAD[c(1:60483),]
row.names (LUAD)<-LUADSEnsemb1_ID
LUAD <-LUAD[,-1]
LUAD_lnc <trans_exp(LUAD,lncrna_only=T)
LUAD_m <-trans_exp(LUAD,mrna_only=T)
pd <-fread("TCGA-LUAD.GDCphenotype.tsv.gz")
sur <fread("TCGA-LUAD.survival.tsv")
group <-make_tcga_group(LUAD_m)
ids_m <-as.data.frame(rownames (LUAD_m))
ids_mSv2 <-rownames (LUAD_m)
degs_m <-get_deg_all(
  LUAD_m,
  group,
  ids_m,
  logFC_cutoff 1,
  ovalue cutoff =0.05.