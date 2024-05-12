
#1. 准备工作


if (!requireNamespace("ggpubr", quietly = TRUE))
  install.packages("ggpubr")  # 安装ggpubr包

if (!requireNamespace("ggbeeswarm", quietly = TRUE))
  install.packages("ggbeeswarm")  # 安装ggpubr包
if (!requireNamespace("ggstatsplot", quietly = TRUE))
  install.packages("ggstatsplot")  # 安装ggpubr包
if (!requireNamespace("PMCMRplus", quietly = TRUE))
  install.packages("PMCMRplus")  # 安装ggpubr包


# 加载必要的库
library(ggplot2)
library(dplyr)
library(ggpubr)  # 加载ggpubr包,加p值
library(ggbeeswarm)#峰群散点


# 使用 mutate() 函数将 cut 列转换为因子变量
# Basic plot
diamonds=c(3,rep(4, 2),rep(5, 4),rep(6, 6),rep(7, 5),rep(8, 2),9)
diamonds=as.data.frame(diamonds)
diamonds[,2]=c(1,rep(2, 2),rep(3, 4),rep(4, 6),rep(5, 5),rep(6, 2),7)
diamonds[,3]=c(5,rep(6, 2),rep(7, 4),rep(8, 6),rep(9, 5),rep(10, 2),11)
diamonds[,4]=c(5,rep(6, 2),rep(7, 4),rep(8, 6),rep(9, 5),rep(10, 2),11)
colnames(diamonds)=c("group1","group2","group3","group4")
diamonds$Sample="sample"


# 将数据从宽格式转换为长格式，以便绘制箱线图
library(tidyr)
re_long <- pivot_longer(diamonds, cols = -Sample, names_to = "group", values_to = "price")
diamonds<-re_long
diamonds$group=as.factor(diamonds$group)





# 2.1 作图
#定义文件名和初始索引
filename <- "箱式小提琴"
index <- 1
# 循环直到找到不存在的文件名
while (file.exists(paste0(filename, "-", index, ".pdf"))) {
  index <- index + 1
}
# 创建新的 PDF 文件
pdf(paste0(filename, "-", index, ".pdf"), onefile = FALSE,width =10, height = 7.5)


# 计算每组数据的平均值和数量
mean_values <- diamonds %>%
  group_by(group) %>%
  summarize(mean_price = mean(price))
count_values <- diamonds %>%
  group_by(group) %>%
  summarize(n_count = n())
# 创建基础图形
set.seed(123)
ggplot(diamonds, aes(group, price, fill = group)) + 
  geom_violin(fill=NA,draw_quantiles = c(0, 1),position = position_dodge(width = 0.75), trim = FALSE,colour="black") +  # 设置小提琴图的填充色为无色，并显示顶部和底部的尖尖
  
  #aes用来传递非颜色量（按分组填充），fill参数用于指定填充颜色（），而color/colour参数用于指定边框颜色（在此为"black"）
  #设置aes后可对fill/color/colour设置全局颜色 scale_fill_manual/ scale_colour_manual，全局设置只对aes生效，对单独的fill/color/colour无效
  #要在小提琴图中忽略离群值，可以在创建geom_violin图层时设置trim参数为TRUE。这将会移除离群值，只展示数据的中间部分。
  geom_boxplot(fill=NA,outlier.shape = NA,width = 0.5, position = position_dodge(width = 0.75)) +  # 添加箱线图
  #在geom_boxplot中，点点通常代表离群值（outliers）,可设置为NA
  
  geom_jitter(shape = 21, size = 4,colour = "black", alpha = 0.8 ,position = position_jitterdodge(dodge.width = 0.75, jitter.width = 2)) +  # 数据点在每一组小提琴图的边框内部，形状为圆形，大小为4，颜色为黑色，透明度为0.5，左右平均分布
  scale_fill_manual(values = c("#00C1A3" ,"#EA8331" , "#9590FF",  "#FF62BC"  ))+  # 设置数据点的填充颜色为青色、红色和紫色
  stat_summary(fun.y = "mean",aes(group = group), geom = "point", shape = 21, size = 6, fill = "firebrick") +  # 添加每个组的平均值点，形状为23，大小为4，填充颜色为红色
  
  geom_text(data = mean_values, aes(label = round(mean_price, 2), x = group, y = mean_price), hjust = -1.8, vjust = 0.5, color = "black", size = 6) +  # 标识每组数据的平均值，文字位于平均值上方，颜色为黑色，大小为3.5，水平显示
  geom_segment(data = mean_values, aes(x = as.numeric(group), xend = as.numeric(group) + 0.35, y = mean_price, yend = mean_price), color = "black", size = 0.5,linetype = "dashed") +
  geom_text(data = count_values, aes(label = paste0("n = ", n_count), x = group, y = 0), color = "black", size = 6, hjust = 0.5, vjust = 2.5) +  # 标注每组数据的数量，文字位于 x 轴下方，颜色为黑色，大小为3.5，水平显示
  stat_compare_means(method = "t.test", label = "p.signif", hide.ns = FALSE) +  # 添加p值标注,ref.group,comparisons可以进行编辑
  theme_light() +  # 设置主题为白色背景

  theme(legend.position = "top",
        text = element_text(size = 40),
       
        panel.background = element_rect(fill = "transparent"), 
        legend.background = element_rect(fill = "transparent"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 关闭 PDF 文件
dev.off()





# 2.2 作图
set.seed(123)
mycolor <- c("#FF9999","#99CC00","#FF9900","#FF99CC","#99CC00","#c77cff")
ggplot(diamonds, aes(group, price, fill = group)) + 
  geom_violin(aes(fill= group,colour=group),width= 0.7,alpha=0.7,draw_quantiles = c(0, 1),position = position_dodge(width = 0.5), trim = FALSE) +  # 设置小提琴图的填充色为无色，并显示顶部和底部的尖尖

  #要在小提琴图中忽略离群值，可以在创建geom_violin图层时设置trim参数为TRUE。这将会移除离群值，只展示数据的中间部分。
  geom_boxplot(fill=NA,outlier.shape = NA,width = 0.4, alpha=0.7,position = position_dodge(width = 0.5),color = "grey40") +  # 添加箱线图
  #在geom_boxplot中，点点通常代表离群值（outliers）,可设置为NA
   #geom_beeswarm(color = "grey60",alpha=0.5,size = 3,method = "swarm",cex = 4,show.legend=F)+
  geom_beeswarm(color = "grey60",alpha=0.5,size = 3,method = "center",cex =4,show.legend=F)+
  #geom_quasirandom(aes(colour = group),alpha=0.75,size = 4,width = 0.3,show.legend=F)+
  geom_line(aes(group=Sample,color=group),alpha=4,linewidth=0.2,show.legend=F)+
  
  scale_colour_manual(values=alpha(mycolor,0.9))+
  #fill参数用于指定填充颜色（在此为"turquoise"和"salmon"），而color参数用于指定边框颜色（在此为"black"）
  scale_fill_manual(values=alpha(mycolor,0.9))+
  labs(x="")+
  theme_classic()











  
  
  
#3.1 分组作图

# 更换填充色，设置分面
ggplot(diamonds,aes(group,price,fill=group)) + 
  geom_violin()+ 
  scale_fill_manual(values = c("palevioletred", "peachpuff", "lightgreen","lightblue", "plum")) +
  facet_wrap(.~price3,ncol = 4)#. ~ price3 意味着将数据按照 price3 列的不同取值进行分组，. 表示其他变量不做改变。ncol = 4 表示每行最多显示 4 个子图。设置分面


  
  
  
#3.2 分组作图(升级版)
library(ggstatsplot)
if (require("PMCMRplus")) {#PMCMRplus：用于多重比较检验
  # to get reproducible results from bootstrapping
  set.seed(123)
  library(ggstatsplot)
  library(dplyr, warn.conflicts = FALSE)
  library(ggplot2)
  
  
  
  # the most basic function call
  grouped_ggbetweenstats(data = filter(ggplot2::mpg, drv != "4"), x = year, y = hwy,
                         grouping.var = drv)
  
  # modifying individual plots using `ggplot.component` argument
  grouped_ggbetweenstats(data = filter(movies_long, genre %in% c("Action", "Comedy"),
                                       mpaa %in% c("R", "PG")), x = genre, y = rating, grouping.var = mpaa, ggplot.component = scale_y_continuous(breaks = seq(1,
                                                                                                                                                               9, 1), limits = (c(1, 9))))

  
  grouped_ggbetweenstats(data = diamonds, x = group, y = price,
                         grouping.var = Sample,
                         violin.args=list(fill= NA,width= 0.7,alpha=0.7,draw_quantiles = c(0, 1),position = position_dodge(width = 0.5), trim = TRUE),  # 设置小提琴图的填充色为无色，并显示顶部和底部的尖尖
                         ggtheme=theme_classic(),
                         point.args=list(size = 4,alpha=0.7,position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.75))
                         )
  
  
}


library(scales)
# 获取默认调色板中的颜色
default_colors <- hue_pal()(10)
# 打印颜色的文字描述
default_colors
"#F8766D" "#D89000" "#A3A500" "#39B600" "#00BF7D"   "#9590FF" "#E76BF3" 
