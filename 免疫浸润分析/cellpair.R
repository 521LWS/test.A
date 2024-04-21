

# Install devtools from CRAN.
install.packages("devtools")

if (!requireNamespace("CSkylarL/TimiGP", quietly = TRUE))
  devtools::install_github("CSkylarL/TimiGP")
if (!requireNamespace("data.table", quietly = TRUE))
  BiocManager::install("data.table",dependencies = TRUE)
if (!requireNamespace("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2",dependencies = TRUE)
# Install TimiGP from GitHub.
# Load the package
library(TimiGP)
library(dplyr)
rm(list=ls())



current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# 将工作路径设置为当前代码文件所在的路径
setwd(current_dir)
#1. 选择数据(后续直接运行，net会伴随NET生成在工作目录，network和edge均可导入cytoscape)
pd <- `TCGA-BRCAprocess_TCGA_pd`
pd<-pd[,c(2,1)]
pd<-mutate_all(pd,as.numeric)
rna = `TCGA-BRCAprocessed_counts_m`







#2. Load cell type and marker annotation ----
data("CellType_Bindea2013_cancer")
geneset <- CellType_Bindea2013_cancer
marker <- unique(geneset$Gene)
#3. Preprocess: TimiCheckEvent & TimiPrePropress ---- 临床信息和转录组图谱的预处理

info <- TimiCheckEvent(pd)
rna <- TimiPrePropress(marker = marker,rna = rna,cohort = rownames(info))



#4. Generate marker pair score: TimiGenePair  ----TimiGenePair将捕获任意两个标记对的逻辑关系，并生成标记对得分（MPS）矩阵：
mps <- TimiGenePair(rna)
#5. Perform univariate Cox regression to find the association between marker pair and survival: TimiCOX ----
res <- TimiCOX(mps = mps,info = info,p.adj = "BH")
mps <- res[[1]]
cox_res <- res[[2]]
# This step takes about 20-30 min, the result has been saved in data as examples
Bindea2013c_COX_MP_SKCM06 <- cox_res
save(Bindea2013c_COX_MP_SKCM06, file = "Bindea2013c_COX_MP_SKCM06.rda")
# 6. Generate Directed Gene Network: TimiGeneNetwork  ----有向基因网络
cox_res <- Bindea2013c_COX_MP_SKCM06
# You can use Cytoscape to visualize the network by setting "export = TRUE"
NET <- TimiGeneNetwork(resdata = cox_res,dataset = "Bindea2013_Cancer",export =TRUE, path = "./")




# 7. Generate Cell Interaction Annotation: TimiCellPair ----载入细胞交互注释表
data(CellType_Bindea2013_cancer)
geneset <- CellType_Bindea2013_cancer
cell_pair <- TimiCellPair(geneset = geneset,core = 20)
# 8. Generate background: TimiBG ----
background <- TimiBG(marker.pair = row.names(cox_res))
# 9. Query: Select marker pairs A_B=1 significantly associated with good prognosis ----
GP <- rownames(cox_res)[which(cox_res$QV<0.05)]
# 10. Enrichment Analysis: TimiEnrich ----富集分析
res <- TimiEnrich(gene = GP, background = background, 
                  geneset = cell_pair, p.adj = "BH",core=20)
# Optional: Permutation to generate FDR: TimiPermFFDR ----这里的niter =40, core =1是我的电脑的极限，情况允许的话改为niter =100，core =3
res <- TimiPermFDR(resdata = res, geneset = geneset, gene = GP,
                   background = background, niter =40, core =1)
# This step takes some time. Please be patient
# This has been saved to data as an example
Bindea2013c_enrich <- res
save(Bindea2013c_enrich,file = "Bindea2013c_enrich.rda")






# 11. Visualization: Dot plot of selected cell interaction: TimiDotplot-----
res <- Bindea2013c_enrich
p <- TimiDotplot(resdata = res,select = c(1:10))
print(p)
# 12. Visualization: Chord Diagram of functional interaction: TimiCellChord----
res <- Bindea2013c_enrich

sum(res$Adjust.P.Value < 0.05) # Returns 54
# Cell Chord Diagram(condition = "P.Value",cutoff = 0.05这里根据情况改condition = "Adjust.P.Value",cutoff = 0.05)
TimiCellChord(resdata = res,dataset = "Bindea2013_Cancer",condition = "Adjust.P.Value",cutoff = 0.05)







#12.1 Chord Diagram of marker pairs in select cell interaction
TimiGeneChord(resdata = res,select = 1)



# 13. Generate Directed Cell Network: TimiCellNetwork  ----
# You can use Cytoscape to visualize the network

res <- Bindea2013c_enrich
NET <- TimiCellNetwork(resdata = res,dataset = "Bindea2013_Cancer",export =TRUE, path = "./",condition = "Adjust.P.Value",cutoff = 0.05)



# 14. Calculate favorability score: TimiFS ----
# Visualization: TimiFSBar 
res <- Bindea2013c_enrich
# Calculate
score <- TimiFS(res)
head(score)
p <- TimiFSBar(score,select = c(1:5,(nrow(score)-2):nrow(score)))
p





