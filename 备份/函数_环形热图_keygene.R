# 5. 环形热图(利用df，geo_exp)(核心要点：y轴的值必须是同一向量,顺序可以不一致，但类别必须一致)







prepare<- function() {
  project=""
  EXP=""
  DEG=""
  hu_man=""
  human_db=""
  
  
}  





heatmap_KEGG<- function(project,EXP,DEG,hu_man, human_db,keygene,df) {
  if (missing(hu_man)) 
    hu_man <- "human"
  if (missing(df)) 
    df <- DEG
  
  if (missing(human_db)) 
    human_db <- "org.Hs.eg.db"
  if (!requireNamespace("tidyverse", quietly = TRUE))
    BiocManager::install("tidyverse",dependencies = TRUE)  
  if (!requireNamespace("ggtree", quietly = TRUE))
    BiocManager::install("ggtree",dependencies = TRUE)
  if (!requireNamespace("treeio", quietly = TRUE))
    BiocManager::install("treeio",dependencies = TRUE)
  if (!requireNamespace("ape", quietly = TRUE))
    BiocManager::install("ape",dependencies = TRUE)
  if (!requireNamespace("ggnewscale", quietly = TRUE))
    BiocManager::install("ggnewscale",dependencies = TRUE)
  if (!requireNamespace("ggtreeExtra", quietly = TRUE))
    BiocManager::install("ggtreeExtra",dependencies = TRUE)
  if (!requireNamespace("RColorBrewer", quietly = TRUE))
    BiocManager::install("RColorBrewer",dependencies = TRUE)
  library(tidyverse)
  library(ggtree)
  library(treeio)
  library(ape)
  library(ggnewscale)
  library(ggtreeExtra)
  library(RColorBrewer)
  
  
  #调用函数
  source("函数\\函数_KEGG分析提取核心基因.R")#默认为human，其他需要谨慎
  circ2=KEGG_circ(project ,DEG=df,n=10, hu_man=hu_man, human_db=human_db,keygene=keygene)
  select_gene <- unique(circ2$genes)
  
  # 第二层（这里挑选的gene略有些随意，根据情况改变）
  #根据表B的行索引对表A的行进行排序
  print(select_gene)
  
  #将小鼠基因转化为大写
  DEG<-as.data.frame(DEG)
  EXP<-as.data.frame(EXP)
  DEG$Gene<-toupper(DEG$Gene)
  
  EXP<-rownames_to_column(EXP,"SYMBOL")
  EXP$SYMBOL<-toupper(EXP$SYMBOL)
  source("函数\\函数_data_cleaning.R")#去除NA，去除gene列的重复行取平均，取log
  EXP <- clean_data(data=EXP, gene="SYMBOL")
  
  
  
  df1 <- DEG[match(select_gene, DEG$Gene), ]
  
  df1 <- na.omit(df1)
  
  df1 <- df1[order(abs(df1$logFC),decreasing = TRUE), ]
  
  df1 <- df1[1:75, ]
  df1 <- na.omit(df1)
  
  
  
  exp1<-EXP[match(df1$Gene, rownames(EXP)), ]
  
  exp1 <- as.data.frame(exp1)
  
  # 将原始行名替换为新的行名序列
  colnames(exp1) <- paste0("sample", seq_len(ncol(exp1)))
  
  exp1<-tibble::rownames_to_column(exp1, var = "Gene")
  exp2 <-pivot_longer(exp1,-Gene)
  #第一层
  if (!requireNamespace("magrittr", quietly = TRUE))
    install.packages("magrittr",dependencies = TRUE)
  
  # 将数据框转换为行名为基因的矩阵
  gene_matrix <- as.matrix(exp1)
  rownames(gene_matrix) <- exp1$Gene
  
  # 计算距离矩阵
  distance_matrix <- dist(gene_matrix)
  
  # 使用hclust()函数进行聚类
  tree <- hclust(distance_matrix, method = "average")
  
  
  #第三层
  df3<-df1
  df3$pvalue<-"logFC"
  
  #第4层
  df4<-df1
  
  
  #第5层
  group<-circ2[,c(3,5)]
  group$group<-"group"
  group$group<-group$term
  
  
  pdf(paste0(project,"/",project, "-","heatmap_KEGG环形热图2.pdf"), width = 10, height = 8)
  
  p3=ggtree(tree,branch.length = "none", layout = "circular",
            linetype = 1,size = 0.5, ladderize = T)+
    layout_fan(angle =180)+
    theme(plot.margin=margin(0,1,-7,0,"cm"))+
    geom_tiplab(offset=11,show.legend=FALSE,size=1.8,
                color = "black",starstroke = 0)+
    
    
    
    geom_fruit(data=exp2,geom=geom_tile,
               mapping=aes(y=Gene,x=name,fill=value),
               pwidth=0.6,offset=0.02,
               axis.params=list(axis="x",text.angle=-90,text.size=2,hjust=0))+
    scale_fill_gradientn(colours = alpha(rev(RColorBrewer::brewer.pal(11,"RdBu")), alpha = 0.9))+
    
    
    new_scale_fill()+
    geom_fruit(data=df4,geom=geom_point, mapping=aes(y=Gene, x=logFC,color=-log10(pvalue)), size=1.5,
               pwidth=0.3,offset = 0.8,axis.params=list(axis="x"  #添加x轴
               ),
               grid.params=list(
                 vline=T
               )
               
    )+
    scale_color_gradientn(colours = colorRampPalette(c("skyblue","#FC5C7D"))(50))+
    
    
    
    
    new_scale_fill()+
    geom_fruit(data=group,geom=geom_tile,
               mapping=aes(y=genes,x=group,fill=term),color="white",
               pwidth=0.2,offset=-0.1)+
    scale_fill_manual(values = brewer.pal(10,"BrBG")[1:26])+
    
    
    
    
    new_scale_fill()+
    geom_fruit(data=df3,geom=geom_tile,
               mapping=aes(y=Gene,x=pvalue,fill=logFC),color="white",
               pwidth=2,offset=0.05)+
    scale_fill_gradientn(colours = alpha(colorRampPalette(c("red","white", "blue"))(9), alpha = 0.9))
  
  print(p3)#这一步可以防止某些文件错误
  dev.off()
  
  
}
#调用函数
#source("函数\\函数_环形热图_keygene.R")#默认为human，其他需要谨慎
#heatmap_KEGG(project=PROJECT,EXP=geo_exp,DEG=df,keygene)



