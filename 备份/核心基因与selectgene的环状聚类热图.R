

project=PROJECT
EXP=geo_exp
DEG=select_df
hu_man=hu_man
human_db=human_db

df=df
#两种方案，一种df=df，在已分析通路中作图。一种直接分析select_df的KEGG需要调整minGSSize = 3参数



if (missing(hu_man)) 
  hu_man <- "human"

if (missing(human_db)) 
  human_db <- "org.Hs.eg.db"



if (!requireNamespace("ReactomePA", quietly = TRUE))
  BiocManager::install("ReactomePA",dependencies = TRUE)
if (!requireNamespace(human_db, quietly = TRUE)) 
  BiocManager::install(human_db,dependencies = TRUE)
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
library(RColorBrewer) 
if (missing(hu_man)) 
  hu_man <- "human"
if (missing(df)) 
  df <- DEG
if (missing(human_db)) 
  human_db <- "org.Hs.eg.db"
if (!requireNamespace("tidyverse", quietly = TRUE))
  BiocManager::install("tidyverse",dependencies = TRUE)  
if (!requireNamespace("ggtree", quietly = TRUE))
  BiocManager::install("ggtree",dependencies = TRUE)
if (!requireNamespace("treeio", quietly = TRUE))
  BiocManager::install("treeio",dependencies = TRUE)
if (!requireNamespace("ape", quietly = TRUE))
  BiocManager::install("ape",dependencies = TRUE)
if (!requireNamespace("ggnewscale", quietly = TRUE))
  BiocManager::install("ggnewscale",dependencies = TRUE)
if (!requireNamespace("ggtreeExtra", quietly = TRUE))
  BiocManager::install("ggtreeExtra",dependencies = TRUE)
if (!requireNamespace("RColorBrewer", quietly = TRUE))
  BiocManager::install("RColorBrewer",dependencies = TRUE)
library(tidyverse)
library(ggtree)
library(treeio)
library(ape)
library(ggnewscale)
library(ggtreeExtra)
library(RColorBrewer)
#绘制net图
if (!requireNamespace("devtools", quietly = TRUE))
  install_github('devtools')
library(devtools)
if (!requireNamespace("ievaKer/aPEAR", quietly = TRUE))
  install_github('ievaKer/aPEAR')

library('aPEAR')



# 1.数据处理（确定物种，更改OrgDb=org.Hs.eg.db。更改或检查DEG是否符合要求colnames(DEG)[1] <- "SYMBOL"，colnames(DEG)[2] <- "logFC"）
# 重命名第一列
colnames(df)[1] <- "SYMBOL"
colnames(df)[2] <- "logFC"
# 从第一列提取基因名称
genename=as.character(df[, 1])
gene_map=bitr(genename,fromType="SYMBOL",toType="ENTREZID",OrgDb=human_db)
# 通过符号将gene_map与df合并
gene_map <- inner_join(gene_map, df, by = "SYMBOL")
# 移除含有NA值的行
gene_map <- na.omit(gene_map)
# 获取排序后的索引
sorted_index <- order(gene_map[, 3], decreasing = TRUE)
# 根据排序后的索引重新排列数据框
gene_map <- gene_map[sorted_index, ]

# gene_map=bitr(genename,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Rn.eg.db)这个却只能识别首字母大写的基因，
#由于chord <- chord_dat(circ, Go3, Go4)中的circ是全拼大写，为了能匹配，这里改大写才行，在gene_map前面改大写，无法被 gene_map=bitr(genename,fromType="SYMBOL",toType="ENTREZID",OrgDb=human_db)识别
gene_map$SYMBOL <- toupper(gene_map$SYMBOL)
# 提取logFC值，并将行名设为Entrez ID
geneList <- gene_map[, 3]
names(geneList) <- as.character(gene_map[, 2])


#2 进行KEGG通路的GSEA
KEGG_gseresult <- gseKEGG(geneList,organism = hu_man,exponent = 1.5,minGSSize = 3,nPerm = 1000,maxGSSize = 1000, pvalueCutoff =1)
KEGG_gseresult <- setReadable(KEGG_gseresult ,OrgDb=human_db,keyType ='ENTREZID')

KEGG_gseresult@result<- subset(KEGG_gseresult@result,abs(NES) >= 1 & pvalue< 0.05 & qvalue< 0.5)
KEGG_gseresult@result<-KEGG_gseresult@result[order(abs(KEGG_gseresult@result$pvalue),decreasing = FALSE), ]

KEGG_gseresult@result<-KEGG_gseresult@result[order(abs(KEGG_gseresult@result$qvalue),decreasing = FALSE), ]
KEGG_gseresult@result<-KEGG_gseresult@result[order(abs(KEGG_gseresult@result$NES),decreasing = TRUE), ]









#此处挑选通路

Go1<-KEGG_gseresult@result[,c(1,1,2,8,11)]
#Go1<-Go_Reactomeresult@result[1:100,c(1,1,2,8,11)]
colnames(Go1)<-c("category","ID","term","adj_pval","genes")
# 假设您的数据框为 data，genes 列为需要修改的列
Go1$genes <- gsub("/", ",", Go1$genes)
#genes	A data frame with columns for 'ID', 'logFC'
Go2<-gene_map[c(1,3)]#注意的是geneID是否是一样的格式，
colnames(Go2)<-c("ID", "logFC")
#合并获取最重要的数据框
circ2 <-circle_dat(Go1,Go2)
select_gene <- unique(circ2$genes)



#3.  画图
# 第二层（这里挑选的gene略有些随意，根据情况改变）
#根据表B的行索引对表A的行进行排序
print(select_gene)

#将小鼠基因转化为大写
DEG<-as.data.frame(DEG)
EXP<-as.data.frame(EXP)
DEG$Gene<-toupper(DEG$Gene)

EXP<-rownames_to_column(EXP,"SYMBOL")
EXP$SYMBOL<-toupper(EXP$SYMBOL)
source("函数\\函数_data_cleaning.R")#去除NA，去除gene列的重复行取平均，取log
EXP <- clean_data(data=EXP, gene="SYMBOL")


df1 <- DEG[match(select_gene,DEG$Gene ), ]

df1 <- na.omit(df1)

df1 <- df1[order(abs(df1$logFC),decreasing = TRUE), ]

df1 <- df1[1:75, ]
df1 <- na.omit(df1)



exp1<-EXP[match(df1$Gene, rownames(EXP)), ]

exp1 <- as.data.frame(exp1)

# 将原始行名替换为新的行名序列
colnames(exp1) <- paste0("sample", seq_len(ncol(exp1)))

exp1<-tibble::rownames_to_column(exp1, var = "Gene")
exp2 <-pivot_longer(exp1,-Gene)
#第一层
if (!requireNamespace("magrittr", quietly = TRUE))
  install.packages("magrittr",dependencies = TRUE)

# 将数据框转换为行名为基因的矩阵
gene_matrix <- as.matrix(exp1)
rownames(gene_matrix) <- exp1$Gene

# 计算距离矩阵
distance_matrix <- dist(gene_matrix)

# 使用hclust()函数进行聚类
tree <- hclust(distance_matrix, method = "average")


#第三层
df3<-df1
df3$pvalue<-"logFC"

#第4层
df4<-df1


#第5层
group<-circ2[,c(3,5)]
group$group<-"group"
group$group<-group$term

library(scales)
pdf(paste0(project,"/",project, "-","heatmap_KEGG环形热图4.pdf"), width = 10, height = 8)

p3=ggtree(tree,branch.length = "none", layout = "circular",
          linetype = 1,size = 0.5, ladderize = T)+
  layout_fan(angle =180)+
  theme(plot.margin=margin(0,1,-7,0,"cm"))+
  geom_tiplab(offset=3.8,show.legend=FALSE,size=1.8,
              color = "black",starstroke = 0)+
  
  
  
  geom_fruit(data=exp2,geom=geom_tile,
             mapping=aes(y=Gene,x=name,fill=value),
             pwidth=0.45,offset=0.02,
             axis.params=list(axis="x",text.angle=-90,text.size=2,hjust=0))+
  scale_fill_gradientn(colours = alpha(rev(RColorBrewer::brewer.pal(11,"RdBu")), alpha = 0.9))+
  
  
  new_scale_fill()+
  geom_fruit(data=df4,geom=geom_point, mapping=aes(y=Gene, x=logFC,color=-log10(pvalue)), size=1.5,
             pwidth=0.2,offset = 0.33,axis.params=list(axis="x"  #添加x轴
             ),
             grid.params=list(
               vline=T
             )
             
  )+
  scale_color_gradientn(colours = colorRampPalette(c("skyblue","#FC5C7D"))(50))+
  
  
  
  
  new_scale_fill()+
  geom_fruit(data=group,geom=geom_tile,
             mapping=aes(y=genes,x=group,fill=term),color="white",
             pwidth=0.2,offset=-0.02)+
  scale_fill_manual(values = hue_pal(l = 30)(11))+
  
  
  
  
  new_scale_fill()+
  geom_fruit(data=df3,geom=geom_tile,
             mapping=aes(y=Gene,x=pvalue,fill=logFC),color="white",
             pwidth=1.5,offset=0.05)+
  scale_fill_gradientn(colours = alpha(colorRampPalette(c("red","white", "blue"))(9), alpha = 0.9))

print(p3)#这一步可以防止某些文件错误
dev.off()


