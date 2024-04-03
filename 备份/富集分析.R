if (!requireNamespace("GOplot", quietly = TRUE))
 install.packages("GOplot",dependencies = TRUE)

library(GOplot)
library(tidyr)

#1. 准备工作（四个文件：来自于聚类结果需要固定格式的的Go1，Go2，感兴趣的基因Go3和go过程Go4，需要调整这四个文件）

#terms	A data frame with columns for 'category', 'ID', 'term', adjusted p-value ('adj_pval') and 'genes'
Go1<-Go_gseresult@result[1:100,c(1,2,3,9,12)]
#Go1<-KEGG_gseresult@result[1:100,c(1,1,2,8,11)]
#Go1<-Go_Reactomeresult@result[1:100,c(1,1,2,8,11)]


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
Go3<-sorted_df[1:5000, ]
Go4<-as.character(Go1[1:7, "term"]) 


#2.作图
# 如果您想要绘制包含 'BP' 类别的 GO 条形图，请使用以下代码：
GOBar(subset(circ, category == 'BP'))
# 如果您想要绘制所有类别的 GO 条形图，请使用以下代码：
GOBar(circ,display='multiple')
# 绘制 GO 圆图，显示 circ$ID 列的前 3 个条目
GOCircle(circ, nsub =8,zsc.col = c("red","white","blue"))


# 3. 可秀图
# 如果您想要绘制包含前三个标签的 GO 气泡图，请使用以下代码(要调整labels = 1.6)：
GOBubble(circ, labels = 1.2,ID = F)
GOBubble(circ, display = 'multiple')
GOBubble(circ, ID = F, table.legend = F)
??GOBubble
# 使用 chord_dat() 函数创建 chord 图所需的数据
chord <- chord_dat(circ, Go3, Go4)
# 查看 chord 数据的前几行
head(chord )
# 使用 GOChord() 函数绘制 chord 图，并设置参数(如果报错行数不匹配（KEGG中），只需要把Go4调少就行了，另外颜色数量也要相应调整)
GOChord(chord, space = 0.01, gene.order = 'logFC', gene.space = 0.25, gene.size = 3,lfc.col=c(rgb(1, 0, 0, 0.9), rgb(1, 1, 1, 0)), ribbon.col=brewer.pal(7, "Set3"))

# # GOCluster() 函数调用示例 #1:4最稳定，其他的话要自己修改，不能直接用于函数出图
GOCluster(circ, Go4,clust.by="logFC",lfc.col=c(rgb(1, 0, 0, 0.5), rgb(1, 1, 1, 0.5), rgb(0, 0, 1, 0.5)),term.col=alpha(brewer.pal(length(Go4), "Paired"), 0.6),lfc.space=1.9,lfc.width=0.5,term.space=-2.2, term.width=1.5)

GOCluster1(circ,  Go3, Go4,clust.by="logFC",lfc.col=c(rgb(1, 0, 0, 0.5), rgb(1, 1, 1, 0.5), rgb(0, 0, 1, 0.5)),term.col=alpha(brewer.pal(length(Go4), "Paired"), 0.6),lfc.space=2.5,lfc.width=0.5,term.space=-3, term.width=1.5 )


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

  
