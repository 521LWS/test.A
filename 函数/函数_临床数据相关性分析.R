prepare<- function() {
  project="TCGA-LIHC"
  date=`TCGA-LIHCprocessed_counts_m`
  pd=`TCGA-LIHCnumeric_TCGA_pd`
  pd=`TCGA-LIHCprocess_TCGA_pd`
  select_gene=c("PCK1","ODC1","PINK1","AGK","PCK2","PGK1","USP38","RPP38","ATAD3A","JAK2","MAPKAP1")
  df=`TCGA-LIHC148df`
} 

correlation_gene_pd1<- function(project,date,pd,df, select_gene) {

  library(WGCNA)  
  library(dplyr)
  library(pheatmap)

  
  if (missing(select_gene))  {
   
    datExpr=t(date)
    pd<- mutate_all(pd, as.numeric)
    #分析方式1
    # 获取表B的行索引
    row_indices <- match(rownames(datExpr), rownames(pd))
    # 根据表B的行索引对表A的行进行排序
    processed2_pd <-pd[row_indices, ]
    
    # 计算基因特征相关性(若出现NA说明数据有误，比如都是同一值,该计算方式会主动忽略原数据中的无效值，但空格可能例外)
    geneTraitSignificance <- as.data.frame(cor(datExpr, processed2_pd, use = "p"))
    geneTraitSignificance[is.na(geneTraitSignificance)] <- 0
    
    # 计算基因特征相关性的 p 值
    p <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(processed2_pd)))
    
    
    tmp <- matrix(ifelse(p < 0.0001, "****",ifelse(p < 0.001, "***",ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(is.na(p), "", ""))))), nrow = nrow(p))
    tmp  <-ifelse(is.na(tmp), "",tmp)
    p1 <- pheatmap(t(geneTraitSignificance),
                   display_numbers =t(tmp),
                   angle_col =45,
                   color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
                   border_color = "white",
                   cellwidth = 0.02, 
                   cellheight = 20,
                   width = 7, 
                   height=9.1,
                   treeheight_col = 0,
                   treeheight_row = 0)
    
    
    pdf(paste0(project,"/",project, "-","临床相关性.pdf"), width = 15, height = 8)
    
    print(p1)
    
    dev.off()
    # 提取 geneTraitSignificance 大于 0.4 且 p 值小于 0.001 的基因
    select_gene_gender<-rownames(geneTraitSignificance[abs(geneTraitSignificance[1]) > 0.3 & p[1] < 0.001, ])
    select_gene_Age<-rownames(geneTraitSignificance[abs(geneTraitSignificance[2]) > 0.2 & p[2] < 0.001, ])
    select_gene_survival<-rownames(geneTraitSignificance[abs(geneTraitSignificance[3]) > 0.15 & p[3] < 0.001, ])
    select_gene_tumertype<-rownames(geneTraitSignificance[abs(geneTraitSignificance[4]) > 0.3 & p[4] < 0.001, ])
    select_gene_DiseaseStatus<-rownames(geneTraitSignificance[abs(geneTraitSignificance[5]) > 0.3 & p[5] < 0.001, ])
    
    saveRDS(select_gene_gender,file = paste0(project,"/",project,"_select_gene_gender.rds"))
    saveRDS(select_gene_Age,file = paste0(project,"/",project,"_select_gene_Age.rds"))
    saveRDS(select_gene_survival,file = paste0(project,"/",project,"_select_gene_survival.rds"))
    saveRDS(select_gene_tumertype,file = paste0(project,"/",project,"_select_gene_tumertype.rds"))
    saveRDS(select_gene_DiseaseStatus,file = paste0(project,"/",project,"_select_gene_DiseaseStatus.rds"))
    
    source("函数\\函数_韦恩图.R")
    venn5(project,set1=t(df[1]),set2=select_gene_DiseaseStatus,set3=select_gene_Age, set4=select_gene_gender, set5=select_gene_survival)
    
    
    
    
    } else {
      datExpr=as.data.frame(t(date[select_gene,])) 
      pd<- mutate_all(pd, as.numeric)
      #分析方式1
      # 获取表B的行索引
      row_indices <- match(rownames(datExpr), rownames(pd))
      # 根据表B的行索引对表A的行进行排序
      processed2_pd <-pd[row_indices, ]
      
      # 计算基因特征相关性(若出现NA说明数据有误，比如都是同一值)
      geneTraitSignificance <- as.data.frame(cor(datExpr, processed2_pd, use = "p"))
      geneTraitSignificance[is.na(geneTraitSignificance)] <- 0
      
      # 计算基因特征相关性的 p 值
      p <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(processed2_pd)))
      
      
      tmp <- matrix(ifelse(p < 0.0001, "****",ifelse(p < 0.001, "***",ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ifelse(is.na(p), "", ""))))), nrow = nrow(p))
      tmp  <-ifelse(is.na(tmp), "",tmp)
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
      pdf(paste0(project,"/",project, "-","KEYGENE临床相关性.pdf"), width = 15, height = 8)
      
      print(p1)
      
      
      dev.off()
      
      
      
      
      
      
      #计算相关性
      
      library(corrplot)
      
      mt<- cor(datExpr, datExpr, use = "p")
      #自定义颜色
      col <- colorRampPalette(c("#83ba9e", "white", "#eb4601"))
      #求P值
      test <-cor.mtest(datExpr, menthod = "pearson", conf.level = 0.95)
      #画图
      
      pdf(paste0(project,"/",project, "-","KEYGENE相关性.pdf"), width = 15, height = 8)
      corrplot(mt, method = "circle", col = col(50), 
               tl.col = "black", tl.cex =1, tl.srt = 45,tl.pos = "lt",
               p.mat = test$p, type = 'lower',cl.pos="n",
               sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2,
               insig = 'label_sig', pch.col = 'gray20', 
               order = 'AOE',addgrid.col = "black")
      p3<-corrplot(mt, method = "number", type = "upper",col =col(50), diag = T,
               tl.col = "n", tl.cex = 1, tl.pos = "n",order = 'AOE',
               add = T,addgrid.col = "black")
      
      dev.off()
      
      
      
      
      
      
      
      
      
  
       
    }

}
#source("函数\\函数_临床数据相关性分析.R")
#correlation_gene_pd1(project,date=`TCGA-LIHCprocessed_counts_m`,pd=`TCGA-LIHCnumeric_TCGA_pd`,df, select_gene) 







correlation_gene_pd2<- function(project,date,pd,df, select_gene) {
  library(WGCNA)  
  library(dplyr)
  library(pheatmap)
  
  if (missing(select_gene))  {
    
    datExpr=t(date)
    pd$Disease_Status <-  ifelse(grepl("III", pd$Disease_Status, ignore.case = TRUE), 4, ifelse(grepl("II", pd$Disease_Status, ignore.case = TRUE), 3, ifelse(grepl("VI", pd$Disease_Status, ignore.case = TRUE), 7, ifelse(grepl("normal", pd$Disease_Status, ignore.case = TRUE), 1, ifelse(grepl("IV", pd$Disease_Status, ignore.case = TRUE), 5, ifelse(grepl("V", pd$Disease_Status, ignore.case = TRUE), 6, ifelse(grepl("I", pd$Disease_Status, ignore.case = TRUE), 2, ifelse(pd$Disease_Status == ""|pd$Disease_Status == " *"|is.na(pd$Disease_Status), NA, 8))))))))
    pd1=pd[order(pd[,"Disease_Status"], decreasing = FALSE),]
    # 循环遍历所有列
    column_names <- names(pd1)
    
    for (col_name in column_names) {
      
      # 获取列中的所有字符，并去重
      unique_values <- unique(pd1[[col_name]])
      # 使用factor函数将字符转换为不同的数字(空格也被命名了，如果不删除无效值就不会被命名)
      pd1[[col_name]] <- as.numeric(factor(pd1[[col_name]], levels = unique_values))
    }
    
    #分析方式1
    # 获取表B的行索引
    row_indices <- match(rownames(datExpr), rownames(pd1))
    # 根据表B的行索引对表A的行进行排序
    processed2_pd <-pd1[row_indices, ]
    
    # 计算基因特征相关性
    geneTraitSignificance <- as.data.frame(cor(datExpr, processed2_pd, use = "p"))
    # 计算基因特征相关性的 p 值
    p <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(processed2_pd)))
    
    
    tmp <- matrix(ifelse(p < 0.0001, "****",ifelse(p < 0.001, "***",ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))), nrow = nrow(p))
    
    p2 <- pheatmap(t(geneTraitSignificance),
                   display_numbers =t(tmp),
                   angle_col =45,
                   color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
                   border_color = "white",
                   cellwidth = 0.02, 
                   cellheight = 20,
                   width = 7, 
                   height=9.1,
                   treeheight_col = 0,
                   treeheight_row = 0)
    pdf(paste0(project,"/",project, "-","临床相关性.pdf"), width = 15, height = 8)
    
    print(p2)
    
    
    dev.off() 
    
    select_gene_gender<-rownames(geneTraitSignificance[abs(geneTraitSignificance[1]) > 0.4 & p[1] < 0.001, ])
    select_gene_Age<-rownames(geneTraitSignificance[abs(geneTraitSignificance[2]) > 0.3 & p[2] < 0.001, ])
    select_gene_survival<-rownames(geneTraitSignificance[abs(geneTraitSignificance[3]) > 0.2 & p[3] < 0.1, ])
    select_gene_tumertype<-rownames(geneTraitSignificance[abs(geneTraitSignificance[4]) > 0.4 & p[4] < 0.001, ])
    select_gene_DiseaseStatus<-rownames(geneTraitSignificance[abs(geneTraitSignificance[5]) > 0.4 & p[5] < 0.001, ])
    
    saveRDS(select_gene_gender,file = paste0(project,"/",project,"_select_gene_gender.rds"))
    saveRDS(select_gene_Age,file = paste0(project,"/",project,"_select_gene_Age.rds"))
    saveRDS(select_gene_survival,file = paste0(project,"/",project,"_select_gene_survival.rds"))
    saveRDS(select_gene_tumertype,file = paste0(project,"/",project,"_select_gene_tumertype.rds"))
    saveRDS(select_gene_DiseaseStatus,file = paste0(project,"/",project,"_select_gene_DiseaseStatus.rds"))
    
    source("函数\\函数_韦恩图.R")
    venn5(project,set1=t(df[1]),set2=select_gene_DiseaseStatus,set3=select_gene_Age, set4=select_gene_gender, set5=select_gene_survival)
    
    
    
    
  } else {
    datExpr=t(date[select_gene,])
    pd$Disease_Status <-  ifelse(grepl("VI", pd$Disease_Status, ignore.case = TRUE), 7, ifelse(grepl("V", pd$Disease_Status, ignore.case = TRUE), 6, ifelse(grepl("IV", pd$Disease_Status, ignore.case = TRUE), 5, ifelse(grepl("normal", pd$Disease_Status, ignore.case = TRUE), 1, ifelse(grepl("III", pd$Disease_Status, ignore.case = TRUE), 4, ifelse(grepl("II", pd$Disease_Status, ignore.case = TRUE), 3, ifelse(grepl("I", pd$Disease_Status, ignore.case = TRUE), 2, ifelse(pd$Disease_Status == ""|pd$Disease_Status == " *"|is.na(pd$Disease_Status), NA, 8))))))))
    pd1=pd[order(pd[,"Disease_Status"], decreasing = FALSE),]
    # 循环遍历所有列
    column_names <- names(pd1)
    
    for (col_name in column_names) {
      
      # 获取列中的所有字符，并去重
      unique_values <- unique(pd1[[col_name]])
      # 使用factor函数将字符转换为不同的数字(空格也被命名了，如果不删除无效值就不会被命名)
      pd1[[col_name]] <- as.numeric(factor(pd1[[col_name]], levels = unique_values))
    }
    
    #分析方式1
    # 获取表B的行索引
    row_indices <- match(rownames(datExpr), rownames(pd1))
    # 根据表B的行索引对表A的行进行排序
    processed2_pd <-pd1[row_indices, ]
    
    # 计算基因特征相关性
    geneTraitSignificance <- as.data.frame(cor(datExpr, processed2_pd, use = "p"))
    # 计算基因特征相关性的 p 值
    p <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(processed2_pd)))
    
    
    tmp <- matrix(ifelse(p < 0.0001, "****",ifelse(p < 0.001, "***",ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))), nrow = nrow(p))
    
    p2 <- pheatmap(t(geneTraitSignificance),
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
    pdf(paste0(project,"/",project, "-","临床相关性.pdf"), width = 15, height = 8)
    
    print(p2)
    
    
    dev.off() 
    
    
    
    

 
}
}
#source("函数\\函数_临床数据相关性分析.R")
#correlation_gene_pd2(project,date=`TCGA-LIHCprocessed_counts_m`,pd=`TCGA-LIHCprocess_TCGA_pd`,df, select_gene) 