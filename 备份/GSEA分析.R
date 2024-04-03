



if (!requireNamespace("ReactomePA", quietly = TRUE))
  BiocManager::install("ReactomePA",dependencies = TRUE)
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) 
  BiocManager::install("org.Hs.eg.db",dependencies = TRUE)
if (!requireNamespace("clusterProfiler", quietly = TRUE)) 
  BiocManager::install("clusterProfiler",dependencies = TRUE)
if (!requireNamespace("enrichplot", quietly = TRUE)) 
  BiocManager::install("enrichplot",dependencies = TRUE)
if (!requireNamespace("biomaRt", quietly = TRUE)) 
  BiocManager::install("biomaRt",dependencies = TRUE)
if (!requireNamespace("tidyverse", quietly = TRUE)) 
  BiocManager::install("tidyverse",dependencies = TRUE)
#data.table 更适用于处理大型数据集和需要高效率运算的情况，由于矛盾这两个只能liberary一个
#而 dplyr 则更适合一般性的数据操作和数据清洗，尤其是对于需要清晰易读的代码和初学者而言。选择使用哪个取决于您的具体需求和个人偏好。
library(tidyverse) 
library(ReactomePA)
library(enrichplot)
library(biomaRt)
library(org.Hs.eg.db) 
library(clusterProfiler)
library(dplyr)
library(org.Hs.eg.db)



DEG=`TCGA-LUADdf`



# 1.数据处理（确定物种，更改OrgDb=org.Hs.eg.db。更改或检查DEG是否符合要求colnames(DEG)[1] <- "SYMBOL"，colnames(DEG)[2] <- "logFC"）
# 重命名第一列
colnames(DEG)[1] <- "SYMBOL"
colnames(DEG)[2] <- "logFC"
# 从第一列提取基因名称
genename=as.character(DEG[, 1])
gene_map=bitr(genename,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db)
# 通过符号将gene_map与DEG合并
gene_map <- inner_join(gene_map, DEG, by = "SYMBOL")
# 移除含有NA值的行
gene_map <- na.omit(gene_map)
# 获取排序后的索引
sorted_index <- order(gene_map[, 3], decreasing = TRUE)
# 根据排序后的索引重新排列数据框
gene_map <- gene_map[sorted_index, ]
# 提取logFC值，并将行名设为Entrez ID
geneList <- gene_map[, 3]
names(geneList) <- as.character(gene_map[, 2])


# 2. 计算富集结果
#计算exponent权重调节的是峰值的偏移一定程度会影响平滑程度，nPerm重复次数，大一点更平滑，但耗时
# 进行GO项的GSEA
Go_gseresult <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", ont = "all", nPerm = 10000, minGSSize = 10,maxGSSize = 1000, pvalueCutoff = 1)
Go_gseresult<-setReadable(Go_gseresult,OrgDb=org.Hs.eg.db,keyType ='ENTREZID')
#获取核心基因
Go_gseresult_core<-as.data.frame(Go_gseresult)
Go_gseresult_core<-Go_gseresult_core[c("ID","core_enrichment")]
#一个搜索函数
get_second_column <- function(keyword, data_frame) {
  row <- which(data_frame[, 1] == keyword)
  if (length(row) == 0) {
    return(NULL)
  } else {
    core1<-unlist(strsplit(data_frame[row, 2], "/"))
    core1<-as.data.frame(core1)
    return(core1)
  }
}
Go_gseresult_core1<-get_second_column(keyword="GO:0004984", data_frame=Go_gseresult_core)
#弦图准备
Go1=Go_gseresult@result[c(1,2,3,7,12)]
colnames(Go1)<-c("category","ID","term","adj_pval","genes")
Go2=gene_map[c(2,3)]
colnames(go2)<-c("ID", "logFC")
library(GOplot)
circ <-circle_dat(Go1,go2)




# 进行KEGG通路的GSEA
KEGG_gseresult <- gseKEGG(geneList,organism = "human",exponent = 1.5, nPerm = 1000,maxGSSize = 1000, pvalueCutoff =1)
KEGG_gseresult <-setReadable(KEGG_gseresult ,OrgDb=org.Hs.eg.db,keyType ='ENTREZID')

# 进行Reactome通路的GSEA
Go_Reactomeresult <- gsePathway(geneList,organism = "human",exponent = 1.5, nPerm = 10000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 1)
Go_Reactomeresult<-setReadable(Go_Reactomeresult,OrgDb=org.Hs.eg.db,keyType ='ENTREZID')

# 3.绘图（富集分析）
dotplot(Go_gseresult,color = "pvalue", showCategory = 20, font.size = 8,label_format=80)
ridgeplot(Go_gseresult,showCategory = 20,label_format=100)

dotplot(KEGG_gseresult, showCategory = 20, label_format=100)
ridgeplot(KEGG_gseresult,showCategory = 20,label_format=100)

dotplot(Go_Reactomeresult, showCategory = 10, label_format=100)
ridgeplot(Go_Reactomeresult,showCategory = 10,label_format=100)

# 4.ggplot绘图（富集分析，可改颜色）
if (!requireNamespace("enrichplot", quietly = TRUE)) 
  BiocManager::install("enrichplot",dependencies = TRUE)
library(ggplot2)
# 根据排序后的索引重新排列数据框
data <- head(as.data.frame(Go_gseresult), 20)
data$count<-data$setSize
# 创建ggplot2图表(scale_colour_gradient2时 midpoint = 0.00010035需要调整)
library(RColorBrewer)
mypalette <- colorRampPalette(c("blue","white","red" ))(50)

ggplot(data, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore))) +
  geom_point(aes(color =-log10(pvalue), size = count), shape = 16) +
  scale_colour_gradient(low = mypalette[1], high = mypalette[40]) +
  #scale_colour_gradientn(colors = mypalette)+
  #scale_colour_gradient2(low = mypalette[1], mid = adjustcolor("grey", alpha.f = 0.5), high = mypalette[9], midpoint = 0.00010035)+
  #guides(color = guide_legend(override.aes = list(shape = 15, size = 5)))+  # 设置图例中色块的形状和大小
  theme_minimal() + 
  labs(x = "geneRatio", y = "Pathway") +
  # 控制坐标轴线的外观
  theme(axis.line = element_line(color = "black", size = 0.5, linetype = "solid"))+
  # 控制刻度线的外观
  theme(axis.ticks = element_line(color = "black", size = 0.25, linetype = "solid"),
        axis.ticks.length = unit(0.08, "cm"))+  
  # 控制刻度标签的外观
  theme(axis.text = element_text(color = "black", size = 7, angle = 0, hjust = 0.5, vjust = 0.5))+
  # 添加或移除绘图区域的边框
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))







# 改为色块儿图例
mypalette <- colorRampPalette(c("blue","white","red" ))(50)
ggplot(data, aes(x = enrichmentScore, y = reorder(Description, enrichmentScore))) +
  geom_point(aes(color = -log10(pvalue), size = count), shape = 16, stroke = 1, alpha = 0.7) +
 
  scale_size_continuous(range=c(2,6))+
  scale_colour_gradientn(colors = mypalette, breaks = seq(min(-log10(data$pvalue)), max(-log10(data$pvalue)), length.out = 20)) +
  theme_minimal() + # 改为色块儿图例
  labs(x = "geneRatio", y = "Pathway") +
  guides(
    color = guide_colorsteps(title = "-log10pvalue")  # 手动指定图例中的刻度标签数量
  ) +# 将图例改为渐变色块形式
  theme(
    axis.line = element_line(color = "black", size = 0.5, linetype = "solid"),
    axis.ticks = element_line(color = "black", size = 0.25, linetype = "solid"),
    axis.ticks.length = unit(0.08, "cm"),
    axis.text = element_text(color = "black", size = 7, angle = 0, hjust = 0.5, vjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, size = 1)
  )







# 5.绘图（GSEA）

gseaplot2(Go_gseresult, 1:3,  pvalue_table = TRUE, rel_heights = c(1.618, 0.618, 1),ES_geom = "dot")  # 绘制GO富集分析结果

gseaplot2(KEGG_gseresult, 1:3, pvalue_table = TRUE, rel_heights = c(1.618, 0.618, 1),ES_geom = "line")# 绘制KEGG富集分析结果

gseaplot2(Go_Reactomeresult, 1:3, geneSets = core_genes, pvalue_table = TRUE, rel_heights = c(1.618, 0.618, 1),ES_geom = "dot")  # 绘制Reactome富集分析结果


































































#6.添加基因（K1是各个通路的基因，顺序对应）
if (!requireNamespace("ggrepel", quietly = TRUE))
  BiocManager::install("ggrepel",dependencies = TRUE)
library(ggrepel)
gseaScores <- getFromNamespace("gseaScores", "DOSE")
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  if (is.numeric(geneSetID)) geneSetID <- object@result[geneSetID, "ID"]#允许输入的是数字（geneSetID所在的行）
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify = TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore)) / 20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}
print(geneList)
k1 <-gsInfo(KEGG_gseresult,c("hsa04740"))
colnames(k1)[6]<-"FC"
k1<-k1[k1$position==1,]
k1$label <- gene_map$SYMBOL[as.numeric(k1$x)]
