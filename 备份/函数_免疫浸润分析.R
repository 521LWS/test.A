prepare<- function() {
PROJECT='TCGA-LUAD'
select_gene=select_gene[1:20]
expr=`TCGA-LUADprocessed_counts_m`
geneset = rio::import("Marker.xlsx")
pd<-`TCGA-LUADprocess_TCGA_pd`
Group<- tinyarray::make_tcga_group(expr)
}    
    
  immune_cell<- function(PROJECT,expr,geneset,Group,select_gene) {
    
    #1.准备工作
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
    
    if (!requireNamespace("ggpubr", quietly = TRUE))
      install.packages("ggpubr")  # 安装ggpubr包
    
    
    
    #BiocManager::install("GSVA")
    library(ggplot2)
    library(tinyarray)
    library(GSVA)
    library(dplyr)
    library(WGCNA)
    library(pheatmap)
    library(ggpubr)  # 加载ggpubr包
    #????ϸ?????ȼ???
    #2. 数据处理
    # 根据细胞类型分割基因集
    geneset <- split(geneset$Metagene, geneset$`Cell type`)
    # 显示前几行基因集
    lapply(geneset[1:3], head)
    
    
    # 计算每个基因的标准差2：表示按列进行操作。这里 2 表示按列操作，1 表示按行操作。sum：表示要对每列执行的操作，这里是求和操作。
    gene_std <- apply(expr, 1, sd)
    # 筛选出标准差小于阈值的基因
    invariant_genes <- names(gene_std[gene_std < 0.01])
    # 从数据中删除保持不变的基因
    expr <- expr[!rownames(expr) %in% invariant_genes, ]
    
    
    
    
    

    # 使用 gsva 进行基因集变异分析
    #默认情况下，参数 kcdf 设置为 “Gaussian”，适用于对数转换的 microarray、RNA-seq 的 log-CPMs、log-RPKMs 或 log-TPMs。如果输入的表达矩阵是 RNA-seq 的 raw count 数据，应将 kcdf 设置为 “Poisson”。
    library(GSVA)
    require(GSVA) 
    re <- gsva(expr, geneset, method = "ssgsea", mx.diff = FALSE, verbose = TRUE,abs.ranking=TRUE,kcdf="Gaussian")
    
    
    
    
    
    
# 图1.    绘制箱线图数据处理
    # 将转置后的 re 转换为数据框
    re_df <- as.data.frame(t(re))
    # 添加样品名分组作为一列
    re_df$Group <- Group
    
    # 将数据从长格式转换为宽格式，以便绘制箱线图
    library(tidyr)
    library(ggplot2)
    re_long <- pivot_longer(re_df, cols = -Group, names_to = "Cell", values_to = "Value")
   #计算平均值
     mean_values <- re_long %>%
      group_by(Group, Cell) %>%
      summarise(mean_Value = mean(Value))
    #re_long代表数据表，Group, Cell代表分组，Value代表数据值
    pdf(paste0(PROJECT,"/",PROJECT, "-","免疫浸润分析.pdf"), width = 20, height = 8)
    
   P1= ggplot(re_long, aes(x = Cell, y = Value, fill = Group)) +
      
      labs(x = "Cell", y = "Value", fill = "Group") +
      geom_violin( draw_quantiles = c(0, 1),position = position_dodge(width = 0.75), trim = FALSE,colour="black") +  # 设置小提琴图的填充色为无色，并显示顶部和底部的尖尖
      
      #fill参数用于指定填充颜色（在此为"turquoise"和"salmon"），而color参数用于指定边框颜色（在此为"black"）
      scale_fill_manual(values = c("salmon", "turquoise"))+
      #要在小提琴图中忽略离群值，可以在创建geom_violin图层时设置trim参数为TRUE。这将会移除离群值，只展示数据的中间部分。
      geom_boxplot(outlier.shape = NA,width = 0.5, position = position_dodge(width = 0.75)) +  # 添加箱线图
      scale_fill_manual(values = c("salmon", "turquoise"))+
      #在geom_boxplot中，点点通常代表离群值（outliers）,可设置为NA
      geom_jitter(shape = 21, size = 0.2,colour = "black", alpha = 0.8, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +  # 数据点在每一组小提琴图的边框内部，形状为圆形，大小为4，颜色为黑色，透明度为0.5，左右平均分布
      scale_fill_manual(values = c("turquoise", "salmon", "mediumpurple1" )) +  # 设置数据点的填充颜色为青色、红色和紫色
      stat_summary(fun.y = "mean",aes(group = Group),position = position_dodge(width = 0.75), geom = "point", shape = 21, size = 1.5, fill = "firebrick") +  # 添加每个组的平均值点，形状为23，大小为4，填充颜色为红色
      stat_compare_means(method = "t.test", label = "p.signif", hide.ns = FALSE) +  # 添加p值标注,ref.group,comparisons可以进行编辑
      theme_light() +  # 设置主题为白色背景
      theme(legend.position = "top",
            panel.background = element_rect(fill = "transparent"), 
            legend.background = element_rect(fill = "transparent"))  +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(legend.position = "bottom") # 将图例放在底部
   
    print(P1)
    dev.off()
    
    
    
    
    
 # 图2. 计算相关系数矩阵
    # 计算基因特征相关性
    geneTraitSignificance <- as.data.frame(cor(t(re), t(expr[select_gene,]), use = "p"))
    # 计算基因特征相关性的 p 值
    p <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(t(re))))
    
    # 将 p 值按阈值进行分类
    tmp <- matrix(ifelse(p < 0.001, "***",ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))) , nrow = nrow(p))
    
    # 绘制热图
    pdf(paste0(PROJECT,"/",PROJECT, "-","免疫浸润相关性基因分析.pdf"), width = 20, height = 8)
    
    p2 = pheatmap(t(geneTraitSignificance),
                   display_numbers = t(tmp),
                   angle_col = 45,
                   color = colorRampPalette(c("steelblue", "white", "salmon"))(100),
                   border_color = "ivory",
                   cellwidth = 20, 
                   cellheight = 20,
                   width = 7, 
                   height = 9.1,
                   treeheight_col = 0,
                   treeheight_row = 0)
    
    print(p2)
    dev.off()
    
    
    
# 图3.
    pdf(paste0(PROJECT,"/",PROJECT, "-","免疫浸润patient分析.pdf"), width = 20, height = 8)
    
    annotation_col= as.data.frame(Group)
    colnames(annotation_col)[1]<- "Sample"
    rownames(annotation_col)<-colnames(re)
    pheatmap(re,
             show_colnames = T, # 不展示行名
             cluster_rows = F, # 不对行聚类
             cluster_cols = F, # 不对列聚类
             annotation_col = annotation_col, # 加注释
             cellwidth=1,cellheight=5, # 设置单元格的宽度和高度
             
             fontsize=5) # 字体大小
    
    dev.off()
    

    
# 图4.  
    anti_tumor <- c('Activated CD4 T cell', 'Activated CD8 T cell', 'Central memory CD4 T cell', 'Central memory CD8 T cell', 'Effector memeory CD4 T cell', 'Effector memeory CD8 T cell', 'Type 1 T helper cell', 'Type 17 T helper cell', 'Activated dendritic cell', 'CD56bright natural killer cell', 'Natural killer cell', 'Natural killer T cell')
    pro_tumor <- c('Regulatory T cell', 'Type 2 T helper cell', 'CD56dim natural killer cell', 'Immature dendritic cell', 'Macrophage', 'MDSC', 'Neutrophil', 'Plasmacytoid dendritic cell')
    
    anti<- as.data.frame(re[gsub('^ ','',rownames(re))%in%anti_tumor,])
    pro<- as.data.frame(re[gsub('^ ','',rownames(re))%in%pro_tumor,])
    anti_n<- apply(anti,2,sum)
    pro_n<- apply(pro,2,sum)
    patient <-colnames(re)
    library(ggplot2)
    data <- data.frame(anti=anti_n,pro=pro_n,patient=patient)
    anti_pro<- cor.test(anti_n,pro_n,method='pearson')
    
    
    
  
    
    if (missing(pd))  {
      
      p1<- ggplot(data,aes(x = anti, y = pro),color=patient) +
        labs(x="Anti-tumor immunity", y="Pro-tumor suppression") +
        geom_point(aes(color=patient),size=3)+geom_smooth(method='lm')+
        annotate("text", label=paste0('R=',round(anti_pro$estimate,4),'\n','p<0.001'))+
        theme(legend.title = element_text(size = 2, face = "bold"), 
              legend.text = element_text(size = 0.5),
              legend.key.size = unit(0.1, "cm"))
      
      
      p2<- ggplot(data,aes(x = anti, y = pro),color=patient) +
        labs(x="Anti-tumor immunity", y="Pro-tumor suppression") +
        geom_point(aes(color=Group),size=3)+geom_smooth(method='lm')+
        annotate("text", label=paste0('R=',round(anti_pro$estimate,4),'\n','p<0.001'))+
        theme(legend.title = element_text(size = 10, face = "bold"), 
              legend.text = element_text(size = 10),
              legend.key.size = unit(1, "cm"))
      
      
      
      
      
      
      pdf(paste0(PROJECT,"/",PROJECT, "-","免疫浸润cor.pdf"), width = 20, height = 8)
      
      print(p1)
      
      print(p2)
      dev.off()
      
      
      } else {
        
  
        
        #2.1. 创建分组
        # 获取表B的行索引
        row_indices <- match(colnames(expr), rownames(pd))
        # 根据表B的行索引对表A的行进行排序
        processed2_pd <-pd[row_indices, ]
        group1 <- make_tcga_group(expr)
        
        group2 <- ifelse(stringr::str_detect(processed2_pd$Disease_Status , "VI"), "stage VI", ifelse(stringr::str_detect(processed2_pd$Disease_Status , "V"), "stage V", ifelse(stringr::str_detect(processed2_pd$Disease_Status , "IV"), "stage IV", ifelse(stringr::str_detect(processed2_pd$Disease_Status , "III"), "stage III", ifelse(stringr::str_detect(processed2_pd$Disease_Status , "II"), "stage II", ifelse(stringr::str_detect(processed2_pd$Disease_Status , "I"), "stage I", ifelse(stringr::str_detect(processed2_pd$Disease_Status , "normal"), "normal", "na")))))))
        unique(group2)
        
        group <-paste0(group1," ",group2)
        
        
        p1<- ggplot(data,aes(x = anti, y = pro),color=patient) +
          labs(x="Anti-tumor immunity", y="Pro-tumor suppression") +
          geom_point(aes(color=group2),size=3)+geom_smooth(method='lm')+
          annotate("text", label=paste0('R=',round(anti_pro$estimate,4),'\n','p<0.001'))+
          theme(legend.title = element_text(size = 10, face = "bold"), 
                legend.text = element_text(size = 10),
                legend.key.size = unit(1, "cm"))
        
        
        pdf(paste0(PROJECT,"/",PROJECT, "-","免疫浸润cor.pdf"), width = 20, height = 8)
        
        print(p1)
        
        dev.off()
        
        
      }
  
    

 
    
    
  

    
    
    
    
    
    
    
        return()
    
    
    
    
    
    
    }
  
  
#调用函数
#source("函数\\函数_免疫浸润分析.R")
#immune_cell(PROJECT,expr=as.matrix(expr),geneset,Group=as.character(Group),select_gene)


