
prepare<- function() {
  project="TCGA-LIHC"
  set1=t(`TCGA-LIHC148df`[1])
  
  set2=select_gene_DiseaseStatus
  set3=select_gene_Age
  set4=select_gene_gender
  set5=select_gene_survival
  set6=select_gene_tumertype
  set7="" 
} 







venn2<- function(project,set1,set2) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE))
    install.packages("RColorBrewer")
if (!requireNamespace("ReactomePA", quietly = TRUE))
  install.packages("VennDiagram")
  library(RColorBrewer)
  library (VennDiagram) 
p2=venn.diagram(x=list(set1,set2),
             
             scaled = F, # 根据比例显示大小
             
             alpha= 0.6, #透明度
             
             lwd=1,lty=1,col=c("#FEB2B3",'#B2B2FE'), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=c( "#FEB2B3",'#B2B2FE'), # 填充色 配色https://www.58pic.com/
             
             category.names = c("Set1", "Set2") , #标签名
             
             cat.dist = 0.02, # 标签距离圆圈的远近
             
             cat.pos = -180, # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 2, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
          
             filename=paste0(project,"/",project, "-",'两组.png'),# 文件保存
             disable.logging = TRUE,
             
             imagetype="png",  # 类型（tiff png svg）
             
             resolution = 400,  # 分辨率
           
            
)
}
#调用函数
#source("函数\\函数_韦恩图.R")
#venn2(project,set1,set2)


venn3<- function(project,set1,set2,set3) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE))
    install.packages("RColorBrewer")
  if (!requireNamespace("ReactomePA", quietly = TRUE))
    install.packages("VennDiagram")
  library(RColorBrewer)
  library (VennDiagram) 
P3=venn.diagram(x=list(set1,set2,set3),
             
             scaled = F, # 根据比例显示大小
             
             alpha= 0.7, #透明度
             
             lwd=1,lty=1,col=c('#B3EAFF',"#FEB2B3","#C6EBB3"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
             
             label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
             
             cex = 2, # 数字大小
             
             fontface = "bold",  # 字体粗细；加粗bold
             
             fill=c('#B3EAFF',"#FEB2B3","#C6EBB3"), # 填充色 配色https://www.58pic.com/
             
             category.names = c("Set1", "Set2","Set3") , #标签名
             
             cat.dist = 0.02, # 标签距离圆圈的远近
             
             cat.pos = c(-120, -240, -180), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
             
             cat.cex = 2, #标签字体大小
             
             cat.fontface = "bold",  # 标签字体加粗
             
             cat.col='black' ,   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
             
             cat.default.pos = "outer",  # 标签位置, outer内;text 外
             
             output=TRUE,
             
             filename=paste0(project,"/",project, "-",'三组.png'),# 文件保存
             disable.logging = TRUE,
             
             imagetype="png",  # 类型（tiff png svg）
             
             resolution = 400,  # 分辨率
             
             compression = "lzw"# 压缩算法
             
)
}
#调用函数
#source("函数\\函数_韦恩图.R")
#venn3(project,set1,set2,set3)


venn4<- function(project,set1,set2,set3,set4) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE))
    install.packages("RColorBrewer")
  if (!requireNamespace("ReactomePA", quietly = TRUE))
    install.packages("VennDiagram")
  library(RColorBrewer)
  library (VennDiagram) 
  p4=venn.diagram(x=list(set1,set2,set3,set4),
                  
                  scaled = F, # 根据比例显示大小
                  
                  alpha= 0.6, #透明度
                  
                  lwd=1,lty=1,col=brewer.pal(4, "Set2"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
                  
                  label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
                  
                  cex = 2, # 数字大小
                  
                  fontface = "bold",  # 字体粗细；加粗bold
                  
                  fill=brewer.pal(4, "Set2"), # 填充色 配色https://www.58pic.com/
                  
                  category.names = c("Set1", "Set2","Set3","Set4") , #标签名
                  
                  cat.dist = c(0.2, 0.2, 0.1, 0.1), # 标签距离圆圈的远近
                  
                  cat.pos = c(-20, 20, -20, 20), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
                  
                  cat.cex = 2, #标签字体大小
                  
                  cat.fontface = "bold",  # 标签字体加粗
                  
                  cat.col=brewer.pal(4, "Set2"),   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
                  
                  cat.default.pos = "outer",  # 标签位置, outer内;text 外
                  
                  output=TRUE,
                  
                  filename=paste0(project,"/",project, "-",'四组.png'),# 文件保存
                  disable.logging = TRUE,
                  
                  imagetype="png",  # 类型（tiff png svg）
                  
                  resolution = 400,  # 分辨率
                  
                  compression = "lzw"# 压缩算法
                  
  )
  
  
  
  
}
#调用函数
#source("函数\\函数_韦恩图.R")
#venn4(project,set1,set2,set3,set4)


venn5<- function(project,set1,set2,set3,set4,set5) {
  if (!requireNamespace("RColorBrewer", quietly = TRUE))
    install.packages("RColorBrewer")
  if (!requireNamespace("ReactomePA", quietly = TRUE))
    install.packages("VennDiagram")
  
  library (VennDiagram) 
  library(RColorBrewer)
  p5=venn.diagram(x=list(set1,set2,set3,set4,set5),
                  
                  scaled = F, # 根据比例显示大小
                  
                  alpha= 0.6, #透明度
                  
                  lwd=1,lty=1,col=brewer.pal(5, "Set2"), #圆圈线条粗细、形状、颜色；1 实线, 2 虚线, blank无线条
                  
                  label.col ='black' , # 数字颜色abel.col=c('#FFFFCC','#CCFFFF',......)根据不同颜色显示数值颜色
                  
                  cex = 2, # 数字大小
                  
                  fontface = "bold",  # 字体粗细；加粗bold
                  
                  fill=brewer.pal(5, "Set2"), # 填充色 配色https://www.58pic.com/
                  
                  category.names = c("Set1", "Set2","Set3","Set4","Set5") , #标签名
                  
                  cat.dist = c(0.2, 0.2, 0.2, 0.2, 0.2), # 标签距离圆圈的远近
                  
                  cat.pos = c(0, -10, 240, 120, 20), # 标签相对于圆圈的角度cat.pos = c(-10, 10, 135)
                  
                  cat.cex = 2, #标签字体大小
                  
                  cat.fontface = "bold",  # 标签字体加粗
                  
                  cat.col=brewer.pal(5, "Set2"),   #cat.col=c('#FFFFCC','#CCFFFF',.....)根据相应颜色改变标签颜色
                  
                  cat.default.pos = "outer",  # 标签位置, outer内;text 外
                  
                  filename=paste0(project,"/",project, "-",'五组.png'),# 文件保存
                  disable.logging = TRUE,
                  
                  imagetype="png",  # 类型（tiff png svg）
                  
                  resolution = 400,  # 分辨率
                  
                  compression = "lzw"# 压缩算法
                  
  )
  p5<<-p5

} 

#调用函数
#source("函数\\函数_韦恩图.R")
#venn5(project,set1,set2,set3,set4,set5)
