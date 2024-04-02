library(ape)

# 对样本进行聚类，生成 hclust 对象
sampleTree <- hclust(dist(datExpr), method = "average")

# 将 hclust 对象转换为 phylo 对象
sampleTree2 <- as.phylo(sampleTree2)

# 将 phylo 对象转换为 Newick 格式
newick_tree <- write.tree(sampleTree2)

# 打印 Newick 格式的树
print(sampleTree2)