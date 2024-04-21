rm(list = ls())#董一点生信制作,IC50，半数抑制浓度， 指在一定条件下，化合物或药物能够抑制生物过程或活性的浓度，使得生物过程或活性被抑制50%。通常，IC50用于评估药物的抗生物活性。IC50值越低，说明药物越强效，因为它可以在更低的浓度下抑制目标生物分子的活性。
if (!requireNamespace("oncoPredict", quietly = TRUE))
  BiocManager::install("oncoPredict",dependencies = TRUE)
if (!requireNamespace("data.table", quietly = TRUE))
  BiocManager::install("data.table",dependencies = TRUE)
if (!requireNamespace("ggplot2", quietly = TRUE))
  BiocManager::install("ggplot2",dependencies = TRUE)
if (!requireNamespace("gpubr", quietly = TRUE))
  BiocManager::install("gpubr",dependencies = TRUE)
if (!requireNamespace("tinyarray", quietly = TRUE))
  BiocManager::install("tinyarray",dependencies = TRUE)
if (!requireNamespace("ggsci", quietly = TRUE))
  BiocManager::install("ggsci",dependencies = TRUE)
if (!requireNamespace("limma", quietly = TRUE))
  BiocManager::install("limma",dependencies = TRUE)
if (!requireNamespace("stringr", quietly = TRUE))
  BiocManager::install("stringr",dependencies = TRUE)

if (!requireNamespace("dplyr", quietly = TRUE))
  BiocManager::install("dplyr",dependencies = TRUE)
if (!requireNamespace("tibble", quietly = TRUE))
  BiocManager::install("tibble",dependencies = TRUE)



library(oncoPredict)
library(data.table)
library(ggplot2)
library(ggpubr)
library(tinyarray)
library(dplyr)
library(tibble)
library(ggsci)
library(limma)
library(stringr)
library(tidyr)
#????ϸ?????ȼ???
# 获取当前代码文件的路径
current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)

# 将工作路径设置为当前代码文件所在的路径
setwd(current_dir)


# 1. 读取并格式化GDSC2 最新训练文件（药物）
GDSCdrug<-function(){
  
  GDSC2_Res<-GDSC2_Res0[c("CELL_LINE_NAME","DRUG_ID","LN_IC50")]
  library(tidyr)
  re <- pivot_wider(GDSC2_Res, names_from = DRUG_ID, values_from = LN_IC50)
  ids<-Cell_Lines_Details_1_[c("COSMIC identifier","Sample Name")]
  ids$`COSMIC identifier`<-paste0("COSMIC_",ids$`COSMIC identifier`)
  names(ids)<-c("COSMIC identifier","CELL_LINE_NAME")
  GDSC2_Res <-as.data.frame(inner_join(re,ids, by = "CELL_LINE_NAME")) 
  rownames(GDSC2_Res)<-GDSC2_Res$`COSMIC identifier`
  GDSC2_Res<-GDSC2_Res[,-1]
  GDSC2_Res<- mutate_all(GDSC2_Res, as.numeric)
  GDSC2_Res <- exp(GDSC2_Res)  #下载的数据是被log转换过的，逆转回去(IC50和LN_IC50的区别)
  A<-as.data.frame(compounds_rel_8_5[c("DRUG_NAME","DRUG_ID")])
  A$DRUG_ID <- as.character(A$DRUG_ID)
  B<-as.data.frame(t(GDSC2_Res)) 
  B<-rownames_to_column(B,"DRUG_ID")
  GDSC2_Res<- inner_join(A,B, by = "DRUG_ID")
  GDSC2_Res <-GDSC2_Res[!duplicated(GDSC2_Res$DRUG_NAME),]
  rownames(GDSC2_Res) <- GDSC2_Res$DRUG_NAME
  
  GDSC2_Res<-t(GDSC2_Res[c(-1,-2)])
  saveRDS(GDSC2_Res,file = "GDSC2_Res.rds")
}











#TCGA此处起一镜到底
#2.准备工作（geo数据需要注意pd）

GDSC2_Expr<-`GDSC2_Expr (RMA Normalized and Log Transformed)`
GDSC2_Res<-GDSC2_Res

exp<-`GSE_exp(mean)`


group1<-c("control","control","control","shYank2","shYank2","shYank2")
processed2_pd<-as.data.frame(group1)
rownames(processed2_pd)<-colnames(exp)
# 设置对比矩阵,control在前
formula <- "shYank2-control"











#geo数据此处以及2.1不运行
pd<-`TCGA-LUADprocess_TCGA_pd`

#2.1. 创建分组
# 获取表B的行索引
row_indices <- match(colnames(exp), rownames(pd))
# 根据表B的行索引对表A的行进行排序
processed2_pd <-pd[row_indices, ]
group1 <- make_tcga_group(exp)

group2 <- ifelse(stringr::str_detect(processed2_pd$Disease_Status , "VI"), "stage VI", ifelse(stringr::str_detect(processed2_pd$Disease_Status , "IV"), "stage IV", ifelse(stringr::str_detect(processed2_pd$Disease_Status , "V"), "stage V", ifelse(stringr::str_detect(processed2_pd$Disease_Status , "III"), "stage III", ifelse(stringr::str_detect(processed2_pd$Disease_Status , "II"), "stage II", ifelse(stringr::str_detect(processed2_pd$Disease_Status , "I"), "stage I", ifelse(stringr::str_detect(processed2_pd$Disease_Status , "normal"), "normal", NA)))))))
unique(group2)

group <-paste0(group1," ",group2)





















#2.2. 处理训练集
# 获取表B的行索引(根据B中从上到下的数据，在A中寻找并获取在A中的位置，组成数字串索引)
row_indices <- match(colnames(GDSC2_Expr), rownames(GDSC2_Res))
# 根据表B的行索引对表A的行进行排序
GDSC2_Res <- GDSC2_Res[row_indices, ]

#删除整行都是缺失值的行
GDSC2_Res <- as.matrix(GDSC2_Res[rowSums(is.na(GDSC2_Res)) != ncol(GDSC2_Res), ]) 

# 获取表B的行索引(根据B中从上到下的数据，在A中寻找并获取在A中的位置，组成数字串索引)
col_indices <- match( rownames(GDSC2_Res),colnames(GDSC2_Expr))
# 根据表B的行索引对表A的行进行排序
GDSC2_Expr<- GDSC2_Expr[, col_indices]
#删除整行都是缺失值的行
GDSC2_Expr <- GDSC2_Expr[rowSums(is.na(GDSC2_Expr)) != ncol(GDSC2_Expr), ]
dim(GDSC2_Res)

identical(rownames(GDSC2_Res),colnames(GDSC2_Expr))


#g10 <- row.names(exp)[1:10]
#exp <- exp[rownames(exp)%in%g10,]

#3.预测
calcPhenotype(trainingExprData = GDSC2_Expr,
              trainingPtype = GDSC2_Res,
              testExprData=as.matrix(exp) ,
              batchCorrect = 'eb',  #   "eb" for ComBat  
              powerTransformPhenotype = TRUE,
              removeLowVaryingGenes = 0.2,
              minNumSamples = 10, 
              printOutput = TRUE, 
              
              removeLowVaringGenesFrom = 'rawData' )
res <- read.csv("calcPhenotype_Output\\DrugPredictions.csv")
#药物名转换前用下边这个
#res <- read.csv("calcPhenotype_Output\\DrugPredictions.csv",header = T , stringsAsFactors = F ,check.names = F)
res <- as.data.frame(res)
rownames(res) <- res$X
res <- res[,-1]



#5..差异分析（探针处理方式为取平均值）(需要更改分组信息)
# 加载所需的包
library(limma)
library(stringr)
res=t(res)
res <- log2(res + 1)
res1=rownames_to_column(as.data.frame(res) )

# 创建设计矩阵
design <- model.matrix(~0 + factor(group1))
colnames(design) <- levels(factor(group1))
print(group1)

# 设置对比矩阵,


contrast.matrix <- makeContrasts(contrasts= formula, levels = design)
print(contrast.matrix)
# 拟合线性模型
fit <- lmFit(res, design)
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
colnames(DEG1)<-c("Drug","logFC","pvalue")
df<- subset(DEG1,logFC <= -0.7& pvalue< 0.05)
df<- df[order(df$logFC,decreasing=FALSE), ]
# 获取表B的行索引(根据B中从上到下的数据，在A中寻找并获取在A中的位置，组成数字串索引)
col_indices <- match(df$Drug,rownames(res))
# 根据表B的行索引对表A的行进行排序
res1<- as.data.frame(t(res[col_indices,])) 

res2<-rownames_to_column(res1,"samples")
processed2_pd<-rownames_to_column(processed2_pd,"samples")




#6.验证

#6.1group1图
for (i in 1:ncol(res1)) {
  Cam <- as.data.frame(res1[,i])
  colnames(Cam) <- "senstivity"
  Cam$Risk <- group1
  boxplot=ggboxplot(Cam, x="Risk", y="senstivity", fill="Risk",
                    xlab="Risk",
                    ylab=paste0(colnames(res1)[i], " senstivity (IC50)"),
                    legend.title="Risk",
                    palette=c("steelblue", "salmon")
  )+
    stat_compare_means(aes(label = ..p.signif..),comparisons = list(c("low risk","high risk")),method="t.test")+#???Ӽ???
    pdf(file=paste0("calcPhenotype_Output\\",colnames(res1)[i], ".pdf"), width=5, height=4.5)
  print(boxplot)
  dev.off()
}
library(ggsci)

processed2_pd$group1 <-group1
processed2_pd=mutate_all(processed2_pd, as.factor)

p3 = res2 %>% select(1:ncol(res2)) %>% 
  inner_join(processed2_pd, "samples") %>% 
  pivot_longer(2:ncol(res2), names_to = "drugs", values_to = "ic50") %>% 
  ggplot(., aes( group1, ic50)) +
  geom_boxplot(aes(fill = group1),
               notch = FALSE, notchwidth = 0.9 # 设置凹槽及凹槽宽度
  ) +
  scale_fill_manual(values = c( "#00AFBB", "#FC4E07")) + # 自定义配色
  
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")+
  facet_wrap(vars(drugs), scales = "free_y", nrow = 2)+
  stat_compare_means()
pdf(file=paste0("calcPhenotype_Output\\", "group1.pdf"), width=15, height=10)
print(p3)
dev.off()








#geo到此结束



#6.2group2图
library(ggsci)

processed2_pd$Disease_Status <-group2
processed2_pd=mutate_all(processed2_pd, as.factor)


p2 = res2 %>% select(1:ncol(res2)) %>% 
  inner_join(processed2_pd, "samples") %>% 
  pivot_longer(2:ncol(res2), names_to = "drugs", values_to = "ic50") %>% 
  ggplot(., aes(Disease_Status, ic50)) +
  geom_boxplot(aes(fill = Disease_Status),
               notch = TRUE, notchwidth = 0.9 # 设置凹槽及凹槽宽度
  ) +
  # scale_fill_manual(values = c("<35" = "#00AFBB", "35-60" = "#E7B800", ">60" = "#FC4E07")) + # 自定义配色
  scale_fill_jama() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")+
  facet_wrap(vars(drugs), scales = "free_y", nrow = 3)+
  stat_compare_means()
pdf(file=paste0("calcPhenotype_Output\\", "group2.pdf"), width=10, height=9)
print(p2)
dev.off()






#6.3group其他图

processed3_pd<-processed2_pd
processed3_pd$Age<-as.numeric(processed3_pd$Age) 
processed3_pd$Age <- ifelse(processed3_pd$Age >= 45, '>45', 
                            ifelse(processed3_pd$Age >= 30, "30-45", "<30"))
processed3_pd$Age  <- factor(processed3_pd$Age , levels = c('<30','30-45','>45'))
processed3_pd <- na.omit(processed3_pd)



p2 = res2 %>% select(1:ncol(res2)) %>% 
  inner_join(processed3_pd, "samples") %>% 
  pivot_longer(2:ncol(res2), names_to = "drugs", values_to = "ic50") %>% 
  ggplot(., aes(x = drugs, y =  ic50,fill=Disease_Status)) +
  geom_boxplot(aes(fill = Disease_Status),
               notch = TRUE, notchwidth = 0.9 # 设置凹槽及凹槽宽度
  ) +
  # scale_fill_manual(values = c("<35" = "#00AFBB", "35-60" = "#E7B800", ">60" = "#FC4E07")) + # 自定义配色
  scale_fill_jama() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none")+
  facet_wrap(vars(Age), scales = "free_y", nrow = 2)+
  stat_compare_means()
pdf(file=paste0("calcPhenotype_Output\\", "group3.pdf"), width=12, height=10)

print(p2)

dev.off()

