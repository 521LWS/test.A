
#1. 准备工作


if (!requireNamespace("ggpubr", quietly = TRUE))
  install.packages("ggpubr")  # 安装ggpubr包
# 加载必要的库
library(ggplot2)
library(dplyr)
library(ggpubr)  # 加载ggpubr包,加p值
# 使用 mutate() 函数将 cut 列转换为因子变量
# Basic plot
diamonds=统计
diamonds$group=as.factor(diamonds$group)

# 2. 作图
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
ggplot(diamonds, aes(group, price, fill = group)) + 
  geom_violin(fill=NA,draw_quantiles = c(0, 1),position = position_dodge(width = 0.75), trim = FALSE,colour="black") +  # 设置小提琴图的填充色为无色，并显示顶部和底部的尖尖
  
  #fill参数用于指定填充颜色（在此为"turquoise"和"salmon"），而color参数用于指定边框颜色（在此为"black"）
  scale_fill_manual(values = c("salmon", "turquoise"))+
  #要在小提琴图中忽略离群值，可以在创建geom_violin图层时设置trim参数为TRUE。这将会移除离群值，只展示数据的中间部分。
  geom_boxplot(fill=NA,outlier.shape = NA,width = 2, position = position_dodge(width = 0.75)) +  # 添加箱线图
  #在geom_boxplot中，点点通常代表离群值（outliers）,可设置为NA
  geom_jitter(shape = 21, size = 4,colour = "black", alpha = 0.8, position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.1)) +  # 数据点在每一组小提琴图的边框内部，形状为圆形，大小为4，颜色为黑色，透明度为0.5，左右平均分布
  scale_fill_manual(values = c("turquoise", "salmon", "mediumpurple1" ))+  # 设置数据点的填充颜色为青色、红色和紫色
  stat_summary(fun.y = "mean",aes(group = group) geom = "point", shape = 21, size = 6, fill = "firebrick") +  # 添加每个组的平均值点，形状为23，大小为4，填充颜色为红色
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

  
  
  
#3. 分组作图

# 更换填充色，设置分面
ggplot(diamonds,aes(group,price,fill=group)) + 
  geom_violin()+ 
  scale_fill_manual(values = c("palevioletred", "peachpuff", "lightgreen","lightblue", "plum")) +
  facet_wrap(.~price3,ncol = 4)#. ~ price3 意味着将数据按照 price3 列的不同取值进行分组，. 表示其他变量不做改变。ncol = 4 表示每行最多显示 4 个子图。设置分面


  
  
  
  
  