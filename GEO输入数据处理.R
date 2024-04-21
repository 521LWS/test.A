
geo <-GSE35958
gse <-"GSE"
PROJECT<-gse
PROJECT<-"GSE35958"
h=180
Formula <- "tumor-normal"
hu_man <- "human"
human_db <- "org.Hs.eg.db"

# 1.数据处理(取探针平均值)
options(scipen = 999)
library(GEOquery)
library(AnnoProbe)
library("tinyarray")


library(plyr)
find_anno(geo$gpl)
ids<- AnnoProbe::idmap(geo$gpl,destdir = tempdir())
geo_exp <- geo$exp
geo_pd <- geo$pd



# 更改格式并命名
geo_exp <- as.data.frame(geo_exp)
geo_exp <- tibble::rownames_to_column(geo_exp, var = "Gene")
# 合并基因表达数据和探针注释
geo_exp <- merge(ids, geo_exp, by.x = "probe_id",by.y ="Gene")
# 按基因符号计算平均表达量
#综合起来，这个函数的作用是对输入的数据框 x 执行列均值计算，但只考虑第三列到最后一列的列。这样，它将返回一个包含从第三列到最后一列的每列的平均值的向量。
geo_exp <- ddply(geo_exp, .(symbol), function(x) colMeans(x[, 3:ncol(x)]))
# 将第一列转换为行名
rownames(geo_exp)<- geo_exp$symbol

library(tidyr)
geo_exp <- separate_rows(geo_exp, symbol, sep = "; ")
geo_exp <- separate_rows(geo_exp, symbol, sep = "/")
# 改为大写(聚类分析不识别小写)
geo_exp$symbol <- toupper(geo_exp$gene)

# 删除原始的ID列
geo_exp <- geo_exp[, -1]



























#此处暂停
# 1.1 样品筛选(关键在于更改样品筛选中的h值和分组信息)
source("函数\\函数_data_cleaning.R")#去除NA，去除gene列的重复行取平均，取log
geo_exp <- clean_data(geo_exp, "symbol")

source("函数\\函数_通过聚类进行样品筛选.R")
process_sample_data =process_sample_data(geo_exp,h)
geo_exp<-process_sample_data
geo_pd <-geo_pd[colnames(geo_exp), ]


library(stringr)

# 1.2 分组（更改分组group_list）
group_list <- ifelse(str_detect(geo_pd$title, "ost"), "tumor", "normal")
# group 变量将被转换为一个因子，并且可以在后续的数据分析中使用
# 打印反向排序后的结果
print(group_list)
tumor_index <- which(group_list  == "tumor")
normal_index <- which(group_list  == "normal")
# 重新排列数据框，确保 "tumor" 和 "normal" 样本连续出现
geo_pd <- geo_pd[c(tumor_index, normal_index), ]

source("函数\\函数_标量化geo临床信息.R")
numeric_geo_pd <- process_clinical_data(geo_pd)

# 获取表B的行索引(根据B中从上到下的数据，在A中寻找并获取在A中的位置，组成数字串索引)
col_indices <- match(rownames(geo_pd), colnames(geo_exp))
# 根据表B的行索引对表A的行进行排序
geo_exp <- geo_exp[,col_indices ]



#1.3 此处需更改
group_list <- ifelse(str_detect(geo_pd$title, "ost"), "tumor", "normal")








#此处暂停

#1.4 芯片数据处理(为防止影响整体代码，包装成了function，只需要更改运行内部代码)

     chuli<-function(geo_exp){
       
       
  hu_man="rat"
        human_db="org.Rn.eg.db"
        
        
        gse <-"GSE"
         PROJECT<-"GSE"
         Formula <- "shYANK2-control"
         h=180
       
  

        geo_exp=as.data.frame(gene_expression_xls[,c(2:8)])
        colnames(geo_exp)[1] <- "gene"

        library(tidyr)
        library(stringr)
        
        geo_exp <- separate_rows(geo_exp, gene, sep = "; ")

        # 1.1 样品筛选(关键在于更改样品筛选中的h值)
        source("函数\\函数_data_cleaning.R")#去除NA，去除gene列的重复行取平均，取log
        geo_exp <- clean_data(data=geo_exp, gene="gene")
        source("函数\\函数_通过聚类进行样品筛选.R")
        process_sample_data =process_sample_data(geo_exp,h)
        geo_exp<-process_sample_data


        group_list=ifelse(str_detect(colnames(geo_exp), "M"), "control", "shYANK2")
        print(group_list)

}









# 2.差异分析（探针处理方式为取平均值）(需要更改分组信息)（limma包处理的是log后的数据，此时得到的2^logFC才接近log前数据的倍数）
# 加载所需的包
library(limma)
library(stringr)

# 创建设计矩阵
design <- model.matrix(~0 + factor(group_list))
colnames(design) <- levels(factor(group_list))
print(group_list )

# 设置对比矩阵,control在前
formula <- Formula

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
DEG1<-DEG[c(1,2,5)]
colnames(DEG1)<-c("Gene","logFC","pvalue")
df<- subset(DEG1,abs(logFC) >= 0 & pvalue< 0.05)


#.3.WGCNA分析
source("函数\\函数_WGCNA出图.R")
WGCNA(project=PROJECT, datExpr=geo_exp, numeric_pd=numeric_geo_pd, module="turquoise")





















#（此处可以暂停，module="turquoise"需要手动更改生成的select-gene需要手动筛选模块儿，可以利用该基因组进行df的筛选）
#df <- df[match(select_gene, df$Gene), ]
#df <- na.omit(df)



#.4.1聚类分析（一般认为|NES|>1，NOM p-val<0.05，FDR q-val<0.25）
source("函数\\函数_聚类分析.R")#默认为human，其他需要谨慎
go_KEGG(project=PROJECT,DEG=df,hu_man=hu_man, human_db=human_db)



# 4.2 环形热图(利用df，geo_exp)

#调用函数
source("函数\\函数_环形热图_依据KEGG筛选的gene.R")#默认为human，其他需要谨慎
heatmap_KEGG(project=PROJECT,EXP=geo_exp,DEG=df,hu_man=hu_man, human_db=human_db)

#4.3免疫浸润分析（需要自己选基因）

select_gene=df[1:20,1]

expr=as.matrix(geo_exp)

geneset = rio::import("Marker.xlsx")
Group<-group_list
Group<-as.character(Group)
#source("函数\\函数_免疫浸润分析.R")
#immune_cell(PROJECT,expr=as.matrix(expr),geneset,Group=as.character(Group),select_gene)

#4.4 临床数据相关性分析（可输入select_gene=c("PCK1")）

geo_exp=`GSE_exp(mean)`
numeric_pd=`numeric_pd`
df=GSEdf

source("函数\\函数_临床数据相关性分析.R")
correlation_gene_pd1(project=PROJECT,date=geo_exp,pd=numeric_pd,df) 


#4.5   
source("函数\\函数_韦恩图.R")
venn5(PROJECT,set1=t(df[1]),set2=TCGA-BRCA_select_gene_DiseaseStatus,set3=TCGA-BRCA_select_gene_Age, set4=TCGA-BRCA_select_gene_gender, set5=TCGA-BRCA_select_gene_survival)
venn2(PROJECT,set1=t(df[1]),set2=select_gene)






#.8.作图
# 探针处理方式不同，以上方式取的均值，在这里取的是最大值，区别相当大，按需求取
# 将mRNA表达数据的行名称转换为数据框形式
ids_m <- as.data.frame(rownames(geo_exp))[c(1,1)]
DEG2<-get_deg(
  
  geo_exp,
  group_list,
  ids_m,
  
  logFC_cutoff = 1,
  pvalue_cutoff = 0.05,
  adjust = FALSE,
  entriz = TRUE,
  species = "human"
)
DEG2<-DEG2[c(8,1,4)]
colnames(DEG2)<-c("Gene","logFC","pvalue")

# 使用get_deg_all()函数找到所有的差异表达基因（DEGs）
group_list <- as.factor(group_list)
DEGS  <-  get_deg_all(
  geo_exp,
  group_list,
  ids_m,
  symmetry = TRUE,
  my_genes = NULL,
  show_rownames = TRUE,
  cluster_cols = TRUE,
  color_volcano = c("#f87669", "grey","#2874C5"),
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


#8.1标签火山图(只需要更改DEG1,注意DEG1是筛选前数据，这样图好看)
if (!requireNamespace('EnhancedVolcano', quietly = TRUE))
  BiocManager::install('EnhancedVolcano')
 library("EnhancedVolcano")
 library(clusterProfiler)
 library(org.Hs.eg.db) 
 
# 重命名第一列

# 整理并重命名
DEG<-DEG[c(1,2,5)]
colnames(DEG)<-c("SYMBOL","logFC","pvalue")
# 从第一列提取基因名称
genename=as.character(DEG[, 1])
gene_map=bitr(genename,fromType="SYMBOL",toType="ENTREZID",OrgDb=human_db)
# 通过符号将gene_map与DEG合并
gene_map <- inner_join(gene_map, DEG, by = "SYMBOL")
# 移除含有NA值的行
gene_map <- na.omit(gene_map)
# 获取排序后的索引
sorted_index <- order(gene_map[, 3], decreasing = TRUE)
# 根据排序后的索引重新排列数据框
gene_map <- gene_map[sorted_index, ]
DEG<-gene_map[c(1,3,4)]
colnames(DEG)<-c("Gene","logFC","pvalue")

p1=EnhancedVolcano(DEG,
                lab = DEG[,1],
                x = 'logFC',
                y = 'pvalue',
             
                xlab = bquote(~Log[2]~ 'fold change'),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                
                #标记设置
               pointSize = c(ifelse(abs(DEG$logFC)>1.5, 8, 2)),  # 使用ifelse语句自定义点的大小
               selectLab =DEG[abs(DEG$logFC)>1.5,1],
               
               
                labSize = 4.0,
                labCol = 'black',
                labFace = 'bold',
                #boxedLabels = TRUE,
                shape = c(1,4,23,25),
               #图例
                colAlpha = 0.8,
                legendLabels=c('Not sig.', 'Log (base 2) FC', 'p-value', 'p-value & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 8,
                legendIconSize = 4.0,
                #箭头
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black',
               
                
                #边框
                gridlines.major = TRUE,
                gridlines.minor = FALSE,
                border = 'full',
                borderWidth = 0.5,
                borderColour = 'black',
                
               
                 )


#此处根据需要修改
p2=p1+
ggplot2::coord_cartesian(xlim=c(-2, 2)) +  # 限制x轴的范围
  
  ggplot2::scale_x_continuous(
    breaks=seq(-2,2, 1))+  # 自定义x轴上的刻度标记
ggplot2::coord_cartesian(ylim=c(0, 8)) +  # 限制y轴的范围
  ggplot2::scale_y_continuous(
    breaks=seq(0,8, 2)) # 自定义x轴上的刻度标记
  
  
  pdf(paste0(gse,"/",gse, "-", "火山图.pdf"), onefile =TRUE,width =9, height =10)

print(p1)

print(p2)

dev.off()
 









#.9.储存

# 创建一个名为 "geo" 的文件夹
dir.create(gse)
# 保存 RDS 文件到新建的文件夹中
filename <- gse
index <- 1
while (file.exists(paste0(filename, "-", index, ".pdf"))) {
  index <- index + 1
}
pdf(paste0(gse,"/",filename, "-", index, ".pdf"), onefile =TRUE,width =15, height =7.50)

DEGS$plots
dev.off()

saveRDS(geo_exp, file = paste0(gse,"/",gse,"_exp(mean).rds"))
saveRDS(DEG1, file = paste0(gse,"/",gse,"DEG1(mean).rds"))
saveRDS(DEG2, file = paste0(gse,"/",gse,"DEG2(max).rds"))
saveRDS(df, file = paste0(gse,"/",gse,"df.rds"))
saveRDS(geo, file = paste0(gse,"/",gse,".rds"))
saveRDS(geo_pd, file = paste0(gse,"/",gse,"_pd.rds"))
saveRDS(numeric_pd, file = paste0(gse,"/",gse," numeric_geo_pd.rds"))
# 存图，需要改名字
#如果执行 pdf() 函数后输出的 PDF 文件是空的，可能是因为在调用 pdf() 函数之后没有进行任何绘图操作，或者绘图操作没有被保存到 PDF 文件中。

