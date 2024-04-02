
if (!requireNamespace("ggtree", quietly = TRUE)) 
  BiocManager::install("ggtree",dependencies = TRUE)
if (!requireNamespace("ggstar", quietly = TRUE)) 
  BiocManager::install("ggstar",dependencies = TRUE)
if (!requireNamespace("ggtreeExtra", quietly = TRUE)) 
  BiocManager::install("ggtreeExtra",dependencies = TRUE)
if (!requireNamespace("treeio", quietly = TRUE)) 
  BiocManager::install("treeio",dependencies = TRUE)
if (!requireNamespace("ggnewscale", quietly = TRUE)) 
  BiocManager::install("ggnewscale",dependencies = TRUE)
if (!requireNamespace("ggsci", quietly = TRUE)) 
  BiocManager::install("ggsci",dependencies = TRUE)
if (!requireNamespace("ape", quietly = TRUE)) 
  BiocManager::install("ape",dependencies = TRUE)


library(ggtree)
library(ggplot2)
library(ggstar)
library(ggtreeExtra)
library(treeio)
library(ggnewscale)
library(ggsci)
library(ape)


#0. 准备工作（打开geo）
source("函数\\函数_geo探针替换.R")
source("函数\\函数_data_cleaning.R")
processed_geo <- clean_data(process_geo(geo))
# 对样本进行聚类(关键是定义sampleTree2)
datExpr <-processed_geo[order(apply(processed_geo, 1, mad), decreasing = TRUE)[1:50], ]
sampleTree <- hclust(dist(datExpr), method = "average")

# 将 hclust 对象转换为 phylo 对象（这个容易看懂，edg代表树枝，第一列与第二列连线）
sampleTree <- as.phylo(sampleTree)
# 使用 ggtree 可视化样本树
par(
  pin = c(12, 9),  # 设置绘图设备的尺寸
  cex = 0.6,       # 设置全局文字缩放因子
  mar = c(0, 4, 2, 0)  # 设置绘图区域的边距
)

#1. 使用 ggtree 可视化样本树
ggtree(sampleTree, layout = "circular") +
  geom_tiplab2(offset = 1, size = 2,vjust = 0.5) +
  geom_text(aes(label = node))  # 您可能想要使用 geom_text 来添加节点标签
  source("函数\\函数_树枝分组函数.R")
  parts <- split_and_fill(df = as.data.frame(as.phylo(sampleTree)$edge), chars = c("52", "65", "53"))


# 如果您想突出显示特定节点并添加标签，您可以尝试以下代码：

#2. 使用 ggtree 可视化样本树
ggtree(sampleTree2, layout = "circular") +
geom_tiplab2(offset = 0.5, size = 2) +
#geom_text(aes(label = node)) + # 如果您想添加节点标签，请取消注释此行
geom_highlight(node = 53, fill = "red") +
geom_highlight(node = 56, fill = "steelblue") +
geom_highlight(node = 72, fill = "green") +
geom_cladelabel(node = 53, label = "virginica", offset = 4, barsize = 25, color = "red", hjust = -0.5,alpha = 0.5) +
geom_cladelabel(node =56, label = "versicolor", offset = 4, barsize = 25, color = "blue", hjust = 1,alpha = 0.5) +
geom_cladelabel(node =72, label = "setosa", offset =4, barsize =25, color = "green", hjust = 3,alpha = 0.5) +
geom_tiplab2(offset = 0.5, size = 2) 

# 3.使用 ggtree 可视化样本树（需要根据树枝分组进行配色，调用树枝分组函数）

ggtree(sampleTree2, layout = "circular") +
  geom_tiplab2(offset = 0.5, size = 2) +
  #geom_text(aes(label = node)) + # 如果您想添加节点标签，请取消注释此行
geom_tree(aes(color = ifelse(node %in% parts[[1]]$V2, "part1",
                             ifelse(node %in% parts[[2]]$V2, "part2",
                                    ifelse(node %in% parts[[3]]$V2, "part3",
                                           ifelse(node %in% parts[[4]]$V2, "part4", "black"))))), size = 1.5) +
  scale_color_manual(values = c("part1" = rgb(0, 0, 1, alpha = 0.6),   # 50% 透明度的蓝色
                                "part2" = "steelblue",
                                "part3" = rgb(0, 1, 0, alpha = 0.6),  # 50% 透明度的绿色
                                "part4" = rgb(1, 0, 0, alpha = 0.6),  # 50% 透明度的红色
                                "other" = rgb(0.5, 0.5, 0.5, alpha = 0.6)))+ # 50% 透明度的灰色

  geom_cladelabel(node = 52, label = "virginica", offset = 4, barsize = 25, color = "blue", hjust =1, alpha = 0.4) +
  geom_cladelabel(node =65, label = "versicolor", offset = 4, barsize = 25, color = "green", hjust = 1,alpha = 0.4) +
  geom_cladelabel(node =53, label = "setosa", offset =4, barsize =25, color = "red", hjust = -1,alpha = 0.4) +
  geom_tiplab2(offset = 0.5, size = 2)+
  theme(legend.position = "none")  # 移除图例

rm(list = ls())  #删除所有对象
