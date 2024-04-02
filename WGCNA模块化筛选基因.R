if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("WGCNA", quietly = TRUE)) 
  BiocManager::install("WGCNA",dependencies = TRUE)
if (!requireNamespace("tinyarray", quietly = TRUE)) 
  BiocManager::install("tinyarray",dependencies = TRUE)
if (!requireNamespace("AnnoProbe", quietly = TRUE)) 
  BiocManager::install("AnnoProbe",dependencies = TRUE)

#1. 准备工作（这里默认载入了geo）
library("WGCNA")
library("tinyarray")
library(AnnoProbe)
source("函数\\函数_标量化geo临床信息.R")
processed_pd <- process_clinical_data(geo$pd)
source("函数\\函数_geo探针替换.R")
processed_geo <- process_geo(geo)
source("函数\\函数_data_cleaning.R")
processed_geo <- clean_data(processed_geo)

## 筛选变化大的基因，并转置
WGCNA_matrix <- processed_geo
## 方法一：SD
##WGCNA_matrix <- t(processed_geo[order(apply(processed_geo, 1, sd), decreasing = TRUE)[1:10000], ])
## 方法二：MAD最常用
WGCNA_matrix <- t(processed_geo[order(apply(processed_geo, 1, mad), decreasing = TRUE)[1:10000], ])
## 方法三：全基因
##WGCNA_matrix <- t(processed_geo)

# 2. 聚类筛选样本(注意临床数据最好一起去除相关样本，设置阅值线和筛选标准cutHeight = 180, minSize = 10)
datExpr <- WGCNA_matrix
sampleNames <- rownames(datExpr)  # 样本名字

sampleTree <- hclust(dist(datExpr), method = "average")  # 对样本进行聚类
par(
  pin = c(12, 9),  # 设置绘图设备的尺寸
  cex = 0.6,       # 设置全局文字缩放因子
  mar = c(0, 4, 2, 0)  # 设置绘图区域的边距
)
xlab = ""       # 设置X轴的标签为空
# 在画布上绘制样本聚类树
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)




# 画一条阅值线，根据实际情况人工设定
abline(h = 180, col = "red")
# 再次对样本进行聚类，确定哪些样本应该被保留,cutHeight参数被设置为160，这意味着将聚类树切成一系列高度为160的子树。
# minSize：minSize参数被设置为10，这意味着每个聚类中至少包含10个样本。
keepSamples <- (cutreeStatic(sampleTree, cutHeight = 180, minSize = 4) ==1)
print(keepSamples)
# 只保留clust=1的样本，即cutreeStatic函数返回的聚类结果中被标记为1的样本
datExpr <- datExpr[keepSamples, ]
# 重新计算样本聚类树
sampleTree <- hclust(dist(datExpr), method = "average")
par(
  pin = c(12, 9),  # 设置绘图设备的尺寸
  cex = 0.6,       # 设置全局文字缩放因子
  mar = c(0, 4, 2, 0)  # 设置绘图区域的边距
)
# 设置 X 轴标签为空
xlab <- ""
plot(sampleTree, main = "Sample clustering after removing outliers", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# 计算基因和样本的数量
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)



# 3. 确定POWER值（只需要修改cex1<-0.85）
# 定义一组用于软阈值的幂次的范围
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# 使用 pickSoftThreshold 函数选择软阈值，verbose = 5 表示输出详细信息
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

cex1<-0.85#不低于0.8
# 将绘图区域分割成 1 行 2 列，便于同时绘制两张图
par(mfrow = c(1, 2))
# 绘制 Scale-Free Topology Model Fit 图
plot(sft$fitIndices[, 1],                # x 轴为软阈值的幂次
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],  # y 轴为 Scale Free Topology Model Fit
     xlab = "Soft Threshold (power)",    # x 轴标签
     ylab = "Scale Free Topology Model Fit, signed RA2",  # y 轴标签
     type = "n",                          # 不绘制点，只绘制图形框架
     main = "Scale independence")         # 图标题
text(sft$fitIndices[, 1],                 # 在图中标注每个点对应的软阈值的幂次
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers,                      # 使用软阈值的幂次作为标签
     cex = cex1,                           # 标签的大小
     col = "red")                          # 标签的颜色
abline(h= cex1, col = "red")              # 在图中添加一条垂直线，表示 Scale-Free 拓扑模型拟合度为 0.85 的阈值
# 绘制 Mean Connectivity 图
plot(sft$fitIndices[, 1],                    # x 轴为软阈值的幂次
     sft$fitIndices[, 5],                    # y 轴为 Mean Connectivity
     xlab = "Soft Threshold (power)",        # x 轴标签
     ylab = "Mean Connectivity",             # y 轴标签
     type = "n",                             # 不绘制点，只绘制图形框架
     main = "Mean Connectivity")             # 图标题
text(sft$fitIndices[, 1],                    # 在图中标注每个点对应的软阈值的幂次
     sft$fitIndices[, 5],
     labels = powers,                        # 使用软阈值的幂次作为标签
     cex = cex1,                             # 标签的大小
     col = "red")                            # 标签的颜色

# 关闭绘图设备
dev.off()

# 输出 powerEstimate 的值
power <- sft$powerEstimate
power



# 4.使用 blockwiseModules 函数进行模块化分析（TOM 矩阵，如果需要，只需要修改power，
#推测：聚类时根据第一步中基因的排序进行聚类，异常表现的基因会被单列出来，只对10000个基因中的5000多个进行了聚类，net$dendrogam[1]储存着树的信息,blockGenes储存着基因的排序编号或者说时索引）
net <- blockwiseModules(
  datExpr,                                # 基因表达数据矩阵
  power = power,                              # 软阈值的幂次
  maxBlockSize = 6000,                    # 最大块大小
  TOMType = "unsigned",                   # TOM 的类型
  minModuleSize = 100,                    # 最小模块大小
  reassignThreshold = 0,                  # 重新指定模块的阈值
  mergeCutHeight = 0.25,                  # 合并模块的高度
  numericLabels = TRUE,                   # 数值标签
  pamRespectsDendro = FALSE,              # PAM 是否尊重树状图
  saveTOMs = FALSE,                        # 是否保存 TOM 矩阵
  verbose = 3                             # 输出详细信息级别
)
# 统计每个模块的基因数目
table(net$colors)
# 将绘图输出保存到文件中，并设置图形大小
png(file = "moduleCluster.png", width = 1200, height = 800)
# 使用 labels2colors 函数将模块颜色转换为标签
mergedColors <- labels2colors(net$colors)
# 使用索引获取参与聚类的基因名称
gene<-rownames(as.data.frame(net$colors))
# 绘制基因树状图和模块颜色
plotDendroAndColors(
  dendro = net$dendrograms[[1]],                  # 基因树状图
  colors = mergedColors[net$blockGenes[[1]]],     # 模块颜色
  groupLabels = "Module colors",                  # 模块颜色标签
  dendroLabels = gene[net$blockGenes[[1]]],                           # 不显示基因树状图标签
  hang = 0.03,                                    # 设置基因树状图的高度
  addGuide = TRUE,                                # 添加颜色图例
  guideHang = 0.05,                               # 设置颜色图例的高度
  main = "Gene dendrogram and module colors"      # 图的主标题
)
# 关闭绘图设备
dev.off()

# 5. 计算模块特征向量与临床特征的相关系数矩阵（需更改临床数据名称,设置图形文字参数）
# 计算模块特征向量（这个特征向量可以理解为一个合成的表达量，能够较好地代表模块内的所有基因的表达情况。）
MEs <-moduleEigengenes(datExpr,mergedColors)$eigengenes
# 按照模块特征向量的顺序进行排序
MEs <- orderMEs(MEs)
# 删除并排列临床数据相应样本
# 获取表B的行索引
row_indices <- match(rownames(datExpr), rownames(processed_pd))
# 根据表B的行索引对表A的行进行排序
processed2_pd <-processed_pd[row_indices, ]






# 计算模块特征向量与临床特征的相关系数矩阵
datTraits <- processed2_pd

moduleTraitCor <- cor(MEs, datTraits, use = "p")
# 计算相关系数的 p 值
nSamples <- ncol(datTraits)
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)
# 将相关系数矩阵和 p 值矩阵写入文件
write.table(file = "step04-modphysiological.cor.xls", moduleTraitCor, sep = "\t", quote = FALSE)
write.table(file = "step04-modphysiological.p.xls", moduleTraitPvalue, sep = "\t", quote = FALSE)
# 绘制模块特征向量与临床特征的热图
pdf(file = "Module_trait_relationships.pdf", width =6, height = 9)
# 创建文本矩阵，包含相关系数和 p 值,\n 是一个转义序列，表示换行符。在字符串中使用 \n 可以实现换行的效果。
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
# 使用 labeledHeatmap 函数绘制热图
par(
  pin = c(6, 9),  # 设置绘图设备的尺寸
  cex = 0.6,       # 设置全局文字缩放因子
  mar = c(4, 8, 4, 4) # 设置绘图区域的边距
)
red_alpha <- adjustcolor("red", alpha.f = 0.5)
blue_alpha <- adjustcolor("blue", alpha.f = 0.5)
white <- "white"
labeledHeatmap(
  Matrix = moduleTraitCor,                      # 相关系数矩阵
  xLabels = colnames(datTraits),                # x 轴标签为临床特征
  yLabels = names(MEs),                         # y 轴标签为模块特征
  ySymbols = names(MEs),                        # y 轴符号为模块特征
  colorLabels = FALSE,                          # 不显示颜色标签
  colors = colorRampPalette(c(red_alpha, white, blue_alpha))(50),  # 设置颜色条
  textMatrix = textMatrix,                      # 文本矩阵
  setStdMargins = TRUE,                         # 设置标准边距
  cex.text = 1.0,                               # 文本大小
  cex.lab = 1.2,                                # 标签大小
  zlim = c(-1, 1),                              # z 轴范围
  main = "Module-trait relationships"           # 图标题
)
# 关闭绘图设备
dev.off()


# 6. 处理数据值（改名任务多）
# 将Disease_stage 转换为数据框
Disease_stage <- as.data.frame(datTraits$Disease_Status)
names(Disease_stage) <- "heart_failure"
# 提取模块名称(只含颜色字符)
modNames <- substring(names(MEs), 3)
# 计算基因模块成员关系
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
# 计算模块成员关系的 p 值
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
# 调整列名
names(geneModuleMembership) <- paste("HF", modNames, sep = "")
names(MMPvalue) <- paste("p.HF", modNames, sep = "")
# 计算基因特征相关性
geneTraitSignificance <- as.data.frame(cor(datExpr, Disease_stage, use = "p"))
# 计算基因特征相关性的 p 值
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
# 调整列名
names(geneTraitSignificance) <- paste("GS.", names(Disease_stage), sep = "")
names(GSPvalue) <- paste("p.GS.", names(Disease_stage), sep = "")


# 7. 选取相关性最强的模块，绘制相关性图
# 选择特定模块（例如棕色模块）
module <- "turquoise"
column <- match(module, modNames)
moduleColors <- mergedColors
moduleGenes <- moduleColors == module

# 保存棕色模块的基因列表
select_gene <- colnames(datExpr)[moduleGenes]
save(select_gene, file = "select_gene.Rdata")

# 绘制模块成员关系 vs. 基因特征关联性的散点图（col = module控制颜色）
pdf(file = "Module_membership_vs_gene_significance.pdf")
par(mfrow = c(1,1))
verboseScatterplot(
  abs(geneModuleMembership[moduleGenes, column]),
  abs(geneTraitSignificance[moduleGenes, 1]),
  xlab = paste("Module Membership in", module, "module"),
  ylab = "Gene significance for Disease",
  main = paste("Module membership vs. gene significance\n"),
  cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module
)
dev.off()




