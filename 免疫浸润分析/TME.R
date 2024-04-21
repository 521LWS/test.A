if (!requireNamespace("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2",dependencies = TRUE)  
if (!requireNamespace("tinyarray", quietly = TRUE))
  BiocManager::install("tinyarray",dependencies = TRUE)
if (!requireNamespace("GSVA", quietly = TRUE))
  BiocManager::install("GSVA",dependencies = TRUE)
if (!requireNamespace("dplyr", quietly = TRUE))
  BiocManager::install("dplyr",dependencies = TRUE)
if (!requireNamespace("WGCNA", quietly = TRUE))
  BiocManager::install("WGCNA",dependencies = TRUE)
if (!requireNamespace("pheatmap", quietly = TRUE))
  BiocManager::install("pheatmap",dependencies = TRUE)


#BiocManager::install("GSVA")
library(ggplot2)
library(tinyarray)
library(GSVA)
library(dplyr)
library(WGCNA)
library(pheatmap)

#????ϸ?????ȼ???
# 获取当前代码文件的路径
current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# 将工作路径设置为当前代码文件所在的路径
setwd(current_dir)




load("choose_gene_1se.Rdata")
load("expr.Rdata")
geneset = rio::import("Marker.xlsx")

geneset = split(geneset$Metagene,geneset$`Cell type`)
lapply(geneset[1:3], head)
re <- gsva(expr, geneset, method="ssgsea",
           mx.diff=FALSE, verbose=FALSE)


draw_boxplot(re,Group,color = c("#1d4a9d","#e5171b"))

#?????????Է???

# 计算基因特征相关性
geneTraitSignificance <- as.data.frame(cor(t(re), t(expr[choose_gene_1se,]), use = "p"))
# 计算基因特征相关性的 p 值
p <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(t(re))))

p[1:4,1:4]

tmp <- matrix(ifelse(p < 0.0001, "****",ifelse(p < 0.001, "***",ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))), nrow = nrow(p))

p1 <- pheatmap(t(geneTraitSignificance),
               display_numbers =t(tmp),
               angle_col =45,
               color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
               border_color = "white",
               cellwidth = 20, 
               cellheight = 20,
               width = 7, 
               height=9.1,
               treeheight_col = 0,
               treeheight_row = 0)


annotation_col= as.data.frame(Group)
colnames(annotation_col)[1]<- "Sample"
rownames(annotation_col)<-colnames(re)
pheatmap(re,
         show_colnames = T, # 不展示行名
         cluster_rows = F, # 不对行聚类
         cluster_cols = TRUE, # 不对列聚类
         annotation_col = annotation_col, # 加注释
         cellwidth=1,cellheight=5, # 设置单元格的宽度和高度
         fontsize=5) # 字体大小



