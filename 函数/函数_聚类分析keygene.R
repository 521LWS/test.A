
prepare<- function() {
  project=""
  DEG=""
  hu_man=""
  human_db=""
  
  
}  





#默认为human，其他需要谨慎
go_KEGG<- function(project,DEG,hu_man, human_db,keygene) {
  if (missing(hu_man))  {
    hu_man <- "human"
  } 
  if (missing(human_db))  {
    human_db <- "org.Hs.eg.db"
  } 
  
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
  
  # 1.数据处理（确定物种，更改OrgDb=org.Hs.eg.db。更改或检查DEG是否符合要求colnames(DEG)[1] <- "SYMBOL"，colnames(DEG)[2] <- "logFC"）
  # 重命名第一列
  colnames(DEG)[1] <- "SYMBOL"
  colnames(DEG)[2] <- "logFC"
  # 从第一列提取基因名称
  genename=as.character(DEG[, 1])
  gene_map=bitr(genename,fromType="SYMBOL",toType="ENTREZID",OrgDb=human_db)
  # 通过符号将gene_map与DEG合并
  gene_map <- inner_join(gene_map, DEG, by = "SYMBOL")
  # 移除含有NA值的行
  gene_map <- na.omit(gene_map)
  # 获取排序后的索引
  sorted_index <- order(gene_map[, 3], decreasing = TRUE)
  # 根据排序后的索引重新排列数据框
  gene_map <- gene_map[sorted_index, ]
  
  # Go_gseresult<- setReadable(Go_gseresult,OrgDb=org.Rn.eg.db,keyType ='ENTREZID')这个转出来的全拼大写基因， gene_map=bitr(genename,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Rn.eg.db)这个却只能识别首字母大写的基因，
  #这里改大写才行，在gene_map前面该大写，无法被 gene_map=bitr(genename,fromType="SYMBOL",toType="ENTREZID",OrgDb=human_db)识别
  gene_map$SYMBOL <- toupper(gene_map$SYMBOL)
  # 提取logFC值，并将行名设为Entrez ID
  geneList <- gene_map[, 3]
  names(geneList) <- as.character(gene_map[, 2])
  
  
  
  
  # 2. 计算富集结果
  #计算exponent权重调节的是峰值的偏移一定程度会影响平滑程度，nPerm重复次数，大一点更平滑，但耗时
  # 2.1 进行GO项的GSEA
  Go_gseresult <- gseGO(geneList, human_db, keyType = "ENTREZID", ont = "all", nPerm = 10000, minGSSize = 10,maxGSSize = 1000, pvalueCutoff = 1)
  Go_gseresult<- setReadable(Go_gseresult,OrgDb=human_db,keyType ='ENTREZID')


  
  
  #弦图准备
  library(GOplot)
  library(tidyr)
  
  #引入改编的函数（1.chord有些特殊，好像是一次性的，这里也要用，加入了函数中。加了标签（角坐标，标签角度问题均已解决））
  GOCluster1<-function (data, Go3,process, metric, clust, clust.by, nlfc, lfc.col, 
                        lfc.space, lfc.width, term.col, term.space,term.width) 
    
  {
    x <- y <- xend <- yend <- width <- space <- logFC <- NULL
    if (missing(metric)) 
      metric <- "euclidean"
    if (missing(clust)) 
      clust <- "average"
    if (missing(clust.by)) 
      clust.by <- "logFC"
    if (missing(nlfc)) 
      nlfc <- 0
    if (missing(lfc.col)) 
      lfc.col <- c(rgb(1, 0, 0, 0.5), rgb(1, 1, 1, 0.5), rgb(0, 0, 1, 0.5))
    # 设置默认值
    if (missing(lfc.space)) {
      lfc.space <- (-0.5)
    } else {
      lfc.space <- lfc.space * (-1)
    }
    
    if (missing(lfc.width)) {
      lfc.width <- (-1.6)
    } else {
      lfc.width <- lfc.space - lfc.width - 0.1
    }
    
    if (missing(term.col)) {
      term.col <- brewer.pal(length(process), "Set3")
    }
    
    if (missing(term.space)) {
      term.space <- lfc.space + lfc.width
    } else {
      term.space <- term.space * (-1) + lfc.width
    }
    
    if (missing(term.width)) {
      term.width <- 2 * lfc.width + term.space
    } else {
      term.width <- term.width * (-1) + term.space
    }
    
    if (clust.by == "logFC") {
      chord <- chord_dat(data, Go3, process)
      distance <- stats::dist(chord[, dim(chord)[2]], method = metric)
    } else {
      chord <- chord_dat(data, Go3, Go4)
      distance <- stats::dist(chord, method = metric)
    }
    
    
    cluster <- stats::hclust(distance, method = clust)
    dendr <- dendro_data(cluster)
    y_range <- range(dendr$segments$y)
    x_pos <- data.frame(x = dendr$label$x, label = as.character(dendr$label$label))
    chord <- as.data.frame(chord)
    chord$label <- as.character(rownames(chord))
    all <- merge(x_pos, chord, by = "label")
    all$label <- as.character(all$label)
    
    if (nlfc) {
      lfc <- all[, c(1, 2, dim(all)[2])]
      lfc$space <- lfc.space
      lfc$width <- lfc.width
    } else {
      lfc <- all[, c(1, 2, dim(all)[2])]
      lfc$space <- lfc.space
      lfc$width <- lfc.width
    }
    
    term <- all[, c(2:(length(process) + 2))]
    color <- NULL
    termx <- NULL
    tspace <- NULL
    twidth <- NULL
    
    for (row in 1:dim(term)[1]) {
      idx <- which(term[row, -1] != 0)
      if (length(idx) != 0) {
        termx <- c(termx, rep(term[row, 1], length(idx)))
        color <- c(color, term.col[idx])
        tmp <- seq(term.space, term.width, length = length(idx) + 1)
        tspace <- c(tspace, tmp[1:(length(tmp) - 1)])
        twidth <- c(twidth, tmp[2:length(tmp)])
      }
    }
    
    # 计算角度
    lfc <- lfc[order(lfc$x), ]
    lfc$angle <- ifelse(seq_along(lfc$x) <= nrow(lfc)/2, -lfc$x * (360/nrow(lfc)) + 90, -lfc$x * (360/nrow(lfc)) - 90)
    
    term_rect <- data.frame(x = termx, width = twidth, space = tspace, col = color)
    legend <- data.frame(x = 1:length(process), label = process)
    ggplot() + geom_segment(data = segment(dendr), aes(x = x,y = y, xend = xend, yend = yend))+
      geom_text(data = lfc, aes(x = x, y =  space+0.5, label = label,angle = angle), size = 2, hjust = 0.5, vjust = 0.5)+  
      geom_rect(data = lfc, aes(xmin = x - 0.5, xmax = x + 0.5, ymin = width, ymax = space, fill = logFC)) + 
      scale_fill_gradient2("logFC", space = "Lab", low = lfc.col[3], mid = lfc.col[2], high = lfc.col[1], 
                           guide = guide_colorbar(title.position = "top", title.hjust = 0.5), 
                           breaks = c(min(lfc$logFC), max(lfc$logFC)), 
                           labels = c(round(min(lfc$logFC)), round(max(lfc$logFC)))) + 
      geom_rect(data = term_rect, aes(xmin = x - 0.5, xmax = x + 0.5, ymin = width, ymax = space), fill = term_rect$col) + 
      geom_point(data = legend, aes(x = x, y = 0.1, size = factor(label, levels = label), shape = NA)) + 
      guides(size = guide_legend("GO Terms", ncol = 3, byrow = TRUE, 
                                 override.aes = list(shape = 22, fill = term.col, size = 8))) + 
      coord_polar() + 
      scale_y_reverse() + 
      
      
      theme_void() + 
      theme(legend.position = "bottom", 
            legend.box = "horizontal", 
            legend.direction = "horizontal")
    
    
  }
  
  
  
  
  
  
  
  Go_gseresult@result<- subset(Go_gseresult@result,abs(NES) >= 1 & pvalue< 0.05 & qvalue< 0.25)
  Go_gseresult@result<-Go_gseresult@result[order(abs(Go_gseresult@result$pvalue),decreasing = FALSE), ]
  Go_gseresult@result<-Go_gseresult@result[order(abs(Go_gseresult@result$qvalue),decreasing = FALSE), ]
  Go_gseresult@result<-Go_gseresult@result[order(abs(Go_gseresult@result$NES),decreasing = TRUE), ]
  #terms	A data frame with columns for 'category', 'ID', 'term', adjusted p-value ('adj_pval') and 'genes'
  Go1<-as.data.frame(Go_gseresult@result[,c(1,2,3,8,12)])

  
  #Go1<-KEGG_gseresult@result[1:100,c(1,1,2,8,11)]
  #Go1<-KEGG_Reactomeresult@result[1:100,c(1,1,2,8,11)]
  colnames(Go1)<-c("category","ID","term","adj_pval","genes")
  # 假设您的数据框为 data，genes 列为需要修改的列
  Go1$genes <- gsub("/", ",", Go1$genes)
  #genes	A data frame with columns for 'ID', 'logFC'
  Go2<-gene_map[c(1,3)]#注意的是geneID是否是一样的格式，
  colnames(Go2)<-c("ID", "logFC")
  #合并获取最重要的数据框
  circ <-circle_dat(Go1,Go2)
  #根据 logFC 的绝对值排序 Go2 数据框的列
  sorted_df <- Go2[rev(order(abs(as.numeric(Go2$logFC)))), ]
  Go3<- subset(sorted_df,abs(logFC) >= 1)
  if (nrow(Go3) < 2000) {
    Go3 <- sorted_df
  } else {
    Go3 <- subset(sorted_df, abs(logFC) >= 1)
  }
  
  

  if (missing(keygene))  {
    if(nrow(Go1)>4){
      n1=4
    } else {
      n1=nrow(Go1)
    }
    Go4<-as.character(Go1[1:n1, 3]) 
  
  } else {
    
    Go_gseresult@result<- Go_gseresult@result[grepl(paste(as.character(keygene), collapse = "|"), Go_gseresult@result[ ,12], ignore.case = TRUE), ]
    Go_gseresult2 <<- Go_gseresult
    if(nrow(Go_gseresult@result)>12){
      n1=12
    } else {
      n1=nrow(Go_gseresult@result)
    }
    Go4<-as.character(Go_gseresult@result[1:n1, 3]) 
 
    }
  print(Go4)
  
  
  # 如果您想要绘制包含前三个标签的 GO 气泡图，请使用以下代码(要调整labels = 1.6)：
  pdf(paste0(project,"/",project, "-","Go_气泡.pdf"), width = 25, height = 8)
  p1<-GOBubble(circ, ID = F, table.legend = F)
  print(p1)
  dev.off()
  
  # 使用 chord_dat() 函数创建 chord 图所需的数据
  pdf(paste0(project,"/",project, "-","Go_弦图2.pdf"), width = 15, height =8)
  chord <- NULL
  chord <- chord_dat(circ, Go3, Go4)
  print(chord)
  par(
    pin = c(12, 9),  # 设置绘图设备的尺寸
    cex = 0.6,       # 设置全局文字缩放因子
    mar = c(1, 1, 1, 1) # 设置绘图区域的边距
  )
  # 使用 GOChord() 函数绘制 chord 图，并设置参数(如果报错选择了未命名的列（KEGG中），只需要把Go4调少就行了，另外颜色数量也要相应调整，这个函数有bug，慎用)
  p1<-GOChord(chord, space = 0.01, gene.order = 'logFC', gene.space = 0.25, gene.size = 3,lfc.col=c(rgb(0, 0, 1, 0.9),rgb(1,1, 1, 0),rgb(1, 0, 0, 0.9)), ribbon.col=brewer.pal(length(Go4), "Set3"))
  print(p1)
  dev.off()
  
  # 可秀图
  pdf(paste0(project,"/",project, "-","Go_聚类图2.pdf"), width =15, height = 8)
  
  par(
    pin = c(12, 9),  # 设置绘图设备的尺寸
    cex = 0.6,       # 设置全局文字缩放因子
    mar = c(1, 1, 1, 1) # 设置绘图区域的边距
  )
  
  pC1<- GOCluster1(circ,  Go3, Go4,clust.by="logFC",lfc.col=c(rgb(0, 0, 1, 0.9),rgb(1,1, 1, 0),rgb(1, 0, 0, 0.9)),term.col=alpha(brewer.pal(length(Go4), "Paired"), 0.6),lfc.space=2.5,lfc.width=0.5,term.space=-3, term.width=1.5 )
  
  print(pC1)
  dev.off()
  
  pdf(paste0(project,"/",project, "-","GSEA-Go_gseresult.pdf"), width = 10, height = 8)
  p1<-gseaplot2(Go_gseresult, 1:n1,  pvalue_table = TRUE, rel_heights = c(1.618, 0.618, 1),ES_geom = "dot")  # 绘制GO富集分析结果
  print(p1)
  dev.off()
  
  # 2. 计算富集结果
  #2. 2 进行KEGG通路的GSEA
  KEGG_gseresult <- gseKEGG(geneList,organism = hu_man,exponent = 1.5, nPerm = 1000,maxGSSize = 1000, pvalueCutoff =1)
  KEGG_gseresult <- setReadable(KEGG_gseresult ,OrgDb=human_db,keyType ='ENTREZID')
 

  KEGG_gseresult@result<-KEGG_gseresult@result[order(abs(KEGG_gseresult@result$pvalue),decreasing = FALSE), ]
  
  KEGG_gseresult@result<-KEGG_gseresult@result[order(abs(KEGG_gseresult@result$qvalue),decreasing = FALSE), ]
  KEGG_gseresult@result<-KEGG_gseresult@result[order(abs(KEGG_gseresult@result$NES),decreasing = TRUE), ]
     Sys.sleep(5)
     
  Go1<-KEGG_gseresult@result[,c(1,1,2,7,11)]
  #Go1<-Go_Reactomeresult@result[1:100,c(1,1,2,7,11)]
  colnames(Go1)<-c("category","ID","term","adj_pval","genes")
  # 假设您的数据框为 data，genes 列为需要修改的列
  Go1$genes <- gsub("/", ",", Go1$genes)
  #genes	A data frame with columns for 'ID', 'logFC'
  Go2<-gene_map[c(1,3)]#注意的是geneID是否是一样的格式，
  colnames(Go2)<-c("ID", "logFC")
  #合并获取最重要的数据框
  circ2 <-circle_dat(Go1,Go2)
  #根据 logFC 的绝对值排序 Go2 数据框的列
  sorted_df <- Go2[rev(order(abs(as.numeric(Go2$logFC)))), ]
  
  Go3<- subset(sorted_df,abs(logFC) >= 1)
  if (nrow(Go3) < 2000) {
    Go3 <- sorted_df
  } else {
    Go3 <- subset(sorted_df, abs(logFC) >= 1)
  }
  
  if (missing(keygene))  {
    if(nrow(Go1)>7){
      n1=7
    } else {
      n1=nrow(Go1)
    }
    Go4<-as.character(Go1[1:n1, 3])  
  } else {
    
    KEGG_gseresult@result<- KEGG_gseresult@result[grepl(paste(as.character(keygene), collapse = "|"), KEGG_gseresult@result[ ,11], ignore.case = TRUE), ]
    KEGG_gseresult2 <<- KEGG_gseresult
    if(nrow(KEGG_gseresult@result)>12){
      n1=12
    } else {
      n1=nrow(KEGG_gseresult@result)
    }
    Go4<-as.character(KEGG_gseresult@result[1:n1, 2]) 
    
  }
print(Go4)
  

  # 可秀图
  chord2 <- NULL
  chord2 <- chord_dat(circ2, Go3, Go4)
  # 使用 chord_dat() 函数创建 chord 图所需的数据
  
  # 继续执行后续的代码
  pdf(paste0(project,"/",project, "-","KEGG_弦图2.pdf"), width = 15, height = 8)
  # 使用 GOChord() 函数绘制 chord 图，并设置参数(如果报错行数不匹配（KEGG中），只需要把Go4调少就行了，另外颜色数量也要相应调整)
  print(chord2)
  par(
    pin = c(12, 9),  # 设置绘图设备的尺寸
    cex = 0.6,       # 设置全局文字缩放因子
    mar = c(1, 1, 1, 1) # 设置绘图区域的边距
  )
  p2<-GOChord(chord2, space = 0.01, gene.order = 'logFC', gene.space = 0.25, gene.size = 3,lfc.col=c(rgb(0, 0, 1, 0.9),rgb(1,1, 1, 0),rgb(1, 0, 0, 0.9)), ribbon.col=brewer.pal(length(Go4), "Set3"))
  print(p2)
  dev.off()
  # 执行某些操作
  pdf(paste0(project,"/",project, "-","KEGG_聚类图2.pdf"), width = 15, height = 8)
  
  par(
    pin = c(12, 9),  # 设置绘图设备的尺寸
    cex = 0.6,       # 设置全局文字缩放因子
    mar = c(1, 1, 1, 1) # 设置绘图区域的边距
  )
  
  pC2<- GOCluster1(data=circ2,  Go3=Go3, process=Go4, clust.by="logFC",lfc.col=c(rgb(0, 0, 1, 0.9),rgb(1,1, 1, 0),rgb(1, 0, 0, 0.9)),term.col=alpha(brewer.pal(length(Go4), "Paired"), 0.6),lfc.space=2.5,lfc.width=0.5,term.space=-3, term.width=1.5 )
  print(pC2)
  dev.off()

  pdf(paste0(project,"/",project, "-","GSEA-Go_KEGG_gseresult.pdf"), width = 10, height = 8)
  p1<-gseaplot2(KEGG_gseresult, 1:n1, pvalue_table = TRUE, rel_heights = c(1.618, 0.618, 1),ES_geom = "line")# 绘制KEGG富集分析结果
  print(p1)
  dev.off()
  # 2. 3 进行Reactome通路的GSEA
  Go_Reactomeresult <- gsePathway(geneList,organism = hu_man,exponent = 1.5, nPerm = 10000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff = 1)
  Go_Reactomeresult <- setReadable(Go_Reactomeresult,OrgDb=human_db,keyType ='ENTREZID')
  
  Go_Reactomeresult@result<- subset(Go_Reactomeresult@result,abs(NES) >= 1 & pvalue< 0.05 & qvalue< 0.25)
  Go_Reactomeresult@result<-Go_Reactomeresult@result[order(abs(Go_Reactomeresult@result$pvalue),decreasing = FALSE), ]
  Go_Reactomeresult@result<-Go_Reactomeresult@result[order(abs(Go_Reactomeresult@result$qvalue),decreasing = FALSE), ]
  Go_Reactomeresult@result<-Go_Reactomeresult@result[order(abs(Go_Reactomeresult@result$NES),decreasing = TRUE), ]
  
  
  Go1<-Go_Reactomeresult@result[,c(1,1,2,7,11)]
  colnames(Go1)<-c("category","ID","term","adj_pval","genes")
  # 假设您的数据框为 data，genes 列为需要修改的列
  Go1$genes <- gsub("/", ",", Go1$genes)
  #genes	A data frame with columns for 'ID', 'logFC'
  Go2<-gene_map[c(1,3)]#注意的是geneID是否是一样的格式，
  colnames(Go2)<-c("ID", "logFC")
  #合并获取最重要的数据框
  circ3 <-circle_dat(Go1,Go2)
  #根据 logFC 的绝对值排序 Go2 数据框的列
  sorted_df <- Go2[rev(order(abs(as.numeric(Go2$logFC)))), ]
  Go3<- subset(sorted_df,abs(logFC) >= 1)
  if (nrow(Go3) < 2000) {
    Go3 <- sorted_df
  } else {
    Go3 <- subset(sorted_df, abs(logFC) >= 1)
  }
  
  
  if (missing(keygene))  {
    if(nrow(Go1)>7){
      n1=7
    } else {
      n1=nrow(Go1)
    }
    Go4<-as.character(Go_Reactomeresult@result[1:n1, 2]) 
  } else {
    Go_Reactomeresult@result<- Go_Reactomeresult@result[grepl(paste(as.character(keygene), collapse = "|"), Go_Reactomeresult@result[ ,11], ignore.case = TRUE), ]
    Go_Reactomeresult2 <<- Go_Reactomeresult
    if(nrow(Go_Reactomeresult@result)>12){
      n1=12
    } else {
      n1=nrow(Go_Reactomeresult@result)
    }
    Go4<-as.character(Go_Reactomeresult@result[1:n1, 2]) 
    
  }
  print(Go4)
  
 
  # 可秀图
  chord3 <- NULL
  chord3 <- chord_dat(circ3, Go3, Go4)
  
  # 执行某些操作
  Sys.sleep(10)  # 暂停5秒钟
  # 继续执行后续的代码
  
  
  
  pdf(paste0(project,"/",project, "-","GSEA-Go_Reactomeresult.pdf"), width = 10, height = 8)
  p1<-gseaplot2(Go_Reactomeresult, 1:n1, pvalue_table = TRUE, rel_heights = c(1.618, 0.618, 1),ES_geom = "line")  # 绘制Reactome富集分析结果
  print(p1)
  dev.off()
  
  
  pdf(paste0(project,"/",project, "-","Go_Reactomeresult_弦图2.pdf"), width = 15, height = 8)
  # 使用 chord_dat() 函数创建 chord 图所需的数据
  par(
    pin = c(12, 9),  # 设置绘图设备的尺寸
    cex = 0.6,       # 设置全局文字缩放因子
    mar = c(1, 1, 1, 1) # 设置绘图区域的边距
  )
  # 使用 GOChord() 函数绘制 chord 图，并设置参数(如果报错行数不匹配（KEGG中），只需要把Go4调少就行了，另外颜色数量也要相应调整)
  
  p3<-GOChord(chord3, space = 0.01, gene.order = 'logFC', gene.space = 0.25, gene.size = 3,lfc.col=c(rgb(0, 0, 1, 0.9),rgb(1,1, 1, 0),rgb(1, 0, 0, 0.9)), ribbon.col=brewer.pal(length(Go4), "Set3"))
  print(p3)
  dev.off()
  
  pdf(paste0(project,"/",project, "-","Go_Reactomeresult_聚类图2.pdf"), width = 15, height = 8)
  
  par(
    pin = c(12, 9),  # 设置绘图设备的尺寸
    cex = 0.6,       # 设置全局文字缩放因子
    mar = c(1, 1, 1, 1) # 设置绘图区域的边距
  )
  
  pC3<- GOCluster1(data=circ3,  Go3=Go3, process=Go4,clust.by="logFC",lfc.col=c(rgb(0, 0, 1, 0.9),rgb(1,1, 1, 0),rgb(1, 0, 0, 0.9)),term.col=alpha(brewer.pal(length(Go4), "Paired"), 0.6),lfc.space=2.5,lfc.width=0.5,term.space=-3, term.width=1.5 )
  
  print(pC3)
  dev.off()
  
  
  # 5.绘图（GSEA）

  
  
  
}
#调用函数
#source("函数\\函数_聚类分析.R")#默认为human，其他需要谨慎
#go_KEGG(project=PROJECT,DEG=df,keygene=C("PCK1","ODC1"))
