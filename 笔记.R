
#1. 中文模式
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")
install.packages("Rcpp",dependencies = TRUE)
install.packages("languageserver")
languageserver::run()

........................................................................
#2. 英文报错
Sys.setenv(LANGUAGE= "en")

# 加入环境变量解决某些包安装问题
system('g++ -v')
system('where make')

.................................................................
#3. 测试相关小技巧
#禁止科学计数法
options(scipen = 999)
#检测Rtool安装与否
system('where make')
#PACKAGE路径
.libPaths()
#程序运行路径

getwd()
setwd()

# 获取当前代码文件的路径
current_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
# 将工作路径设置为当前代码文件所在的路径
setwd(current_dir)

#加载包
library(Rcpp)
#查看占用
search()
#清除占用
detach("package:Rcpp", unload = TRUE)

#清除缓存
rm(ids_570)
rm(list = ls())  # 删除所有对象
remove.packages("Rcpp")

#重新安装
install.packages("languageserver", clean = TRUE, repos=NULL, type="source")

#出现错误时停止后续所有运行并给你选项来试图恢复
options(error = recover)
#可以使用来恢复默认行为
options(error = NULL) 

#获取全部包名称
installed_packages <- installed.packages()
installed_packages[, "Package"]
#重新加载包（可选）：
library()
#或者重新加载已经加载的包：
devtools::reload()

...............................................................
#4 设置github
#设置路径
cd "C:/Users/21989/Documents/test.A"


#设置账户
git config --global user.email "2628220830@qq.com"
git config --global user.name "521LWS"
#设置密码
gpg --full-generate-key #生成密钥
gpg --list-keys #查看密钥
gpg --armor --export 33DBB18B01C92DD677C2AF74AD886BCCCE17DC51 #生成符合github格式的密钥
git config --global user.signingkey 33DBB18B01C92DD677C2AF74AD886BCCCE17DC51
git config  --global credential.helper manager

#处理冲突
git push --set-upstream origin master
git stash
git pull origin main --allow-unrelated-histories
git stash pop
....................................................................
#help
help(package=ggplot2)
??ggplot2
help.search ("volcano")

av <- available.packages()
filtered_pkgs <- av[grep("WGCNA", av[, "Package"]), ]
filtered_pkgs
.....................................................................
....................................................................
#包
#从CRAN
#安装依存
install.packages("WGCNA",dependencies = TRUE)
#从BiocManager#这段代码是在R中用于检查是否已安装BiocManager包的语句。
package_name  <- "BiocManager"
if (!requireNamespace(package_name , quietly = TRUE)) 
  install.packages(package_name )
BiocManager::install("WGCNA",dependencies = TRUE)
#强制安装
BiocManager::install("WGCNA", dependencies = TRUE, force = TRUE)
BiocManager::install("GO.db")
#从github安装,这个一般少用，因为必须知道包来源于谁，比如GO.db来源于Bioconductor
install.packages("devtools")
devtools::install_github("Bioconductor/GO.db")
#当显示没有权限安装和移动包时，最合适的方法就是清除ENVIRONMENT,重启并手动update包后再安装


# 定义要卸载的包名称
package_name <- "WGCNA"
# 找出包及其依赖项
dependencies <- tools::package_dependencies(package_name, recursive = TRUE)$package

# 将包和其依赖项一起卸载
to_remove <- c(package_name, dependencies)
for (pkg in to_remove) {
  if (pkg %in% rownames(installed.packages())) {
    remove.packages(pkg)
  }
}
.....................................................................

#精髓：R是统计学，因此多为向量（集合），而很少标量,c代表连接和位置
x<-c(3:6)
y<-c(2:5)
2*x+y

#变量
a=1:10
dim(a)=c(2,5)
#函数
y<-c(2:5)
seq(from=1,to=100,by=2)
rep(1:5,each=5,times=2)
rep(x,y)
#索引(逻辑索引/位置索引均可,names之后可以用名称索引)
x[x>4 & x>5]
x[-2]#位置
x[c(2,3)]#位置
x %in% c(3)#生成逻辑值
x[x %in% c(3)]#索引逻辑
names(y)<-x
#编辑向量
x[c(2:5)]<-c(2:5)#已改变x
append(x=x,99,3)#未改变x，需赋值
x
............................................................
数据处理分析
#1.移除环境文件
rm(ids, ids_570, geo_pd, geo_exp, fit, fit2, df, data, merged_data,probe_annotation,genes_expr,geo,exp,exp_exp,ensembl,averaged_data,contrast.matrix,DEG,DEGS)
rm(list = ls())  #删除所有对象
rm(net)
# 2..载入
load("TCGA-COAD_RNA.Rdata")
#使用read.csv()函数载入CSV文件
data <- read.csv("your_file.csv", header = TRUE)
# 读取以制表符分隔的文本文件
data <- read.delim("file.txt", header = TRUE)
# 以TXT格式读取基因表达数据
data <- read.table(file_path, header = TRUE, sep = "\t")
# 运行已储存模块
source("data_cleaning.R")



# 4. 查找(切记先进行WGCNA_matrix<-as.data.frame(WGCNA_matrix))
# 筛选
df1<-df[select_gene, ]
# 筛选
protein_coding_rows <- data[data$type == "protein-coding", ]
##转置
WGCNA_matrix <- t(WGCNA_matrix)
#根据关键词提取必要的列
selected_columns <- WGCNA_matrix[, grepl("LEP", names(WGCNA_matrix), ignore.case = TRUE)]
keyword <- c("Patient_ID", "Age", "Sex", "gender","Disease_Status","diagno")
selected_columns <- clinical_data[, grepl(paste(keyword, collapse = "|"), names(clinical_data), ignore.case = TRUE)]
# 将疾病状态转换为数值标量（包含Healthy或者是无效值,ignore.case = TRUE忽略关键词大小写）
selected_columns$Disease_Status <- ifelse(grepl("Healthy", selected_columns$Disease_Status, ignore.case = TRUE) | is.na(selected_columns$Disease_Status), 0, 1)
# 对比两列数据异同
diff_data <- query$submitter_id[!(query$submitter_id %in% clinical$bcr_patient_barcode)]
#找到某一列的某个数值所在的一行
row <- which(DEG[, 1] == "TRAPPC1")
# 使用逻辑索引筛选数据
result <- query[query$submitter_id == "TCGA-5M-AAT5", ]
# 在列名中查找包含"LEP"的列
if ("LEP" %in% colnames(WGCNA_matrix)) {
  print("LEP 列存在于 WGCNA_matrix 中。")
} else {
  print("LEP 列不存在于 WGCNA_matrix 中。")
}
#这将把 genes 列中所有的斜杠替换为逗号。如果您只想修改某些特定行中的值，可以使用子集操作来实现。例如，只想修改前 10 行的值：
data$genes[1:10] <- gsub("/", ",", data$genes[1:10])





# 6.操作行列

gene_expression_data <- read.table("gene_expression_data.txt", header = TRUE)
# 提取列
df<-DEG[,c(1,2,3)]
# 删除列
geo_exp <- geo_exp[, -1]
colnames(df)<-c("Gene","10g2FC","pvalue")
#转置
WGCNA_matrix <- t(fpkm)
# 提取行
processed_geo <- processed_geo[c(1,2,5), ]
df <- DEG[c(1,2,5), 1]
# 行名入表
DEG <- tibble::rownames_to_column(DEG, var = "Gene")
# 将第一列转换为行名，为了方便其他包处理
rownames(geo_exp)<- geo_exp$symbol
# (行列区分就在于逗号在左在右，向量相当于数据表的一行，在公式运用上有一定相似度而不需要指明行列)
#交换第一行和第二行的顺序
tableB <- tableB[c(2, 1), ]

# 获取表B的行索引(根据B中从上到下的数据，在A中寻找并获取在A中的位置，组成数字串索引)
row_indices <- match(rownames(tableB), rownames(tableA))
# 根据表B的行索引对表A的行进行排序
tableA_sorted <- tableA[row_indices, ]
# 根据表B的行索引对表A的行进列排序
processed_geo <- geo_data[, row_indices]
# 将 blockgenes 数据框的第一列的数值作为索引
vecto <- blockgenes[[1]]

# 使用索引从 color 数据框的第一列中提取相应位置的值
dendro$labels <- color[vecto, 1]
a<-DEG[1, 1]

# 7. 处理数据
#重命名
colnames(Go1)<-c("a", "b")
#使用subset函数删除包含空格的行
cleaned_data <- subset(data, !grepl("\\s", rownames(data)))
# 移除含有NA值的行
gene_map <- na.omit(gene_map)
# 假设 data 是包含数据的数据框，且包含一列名为column_name的字符列
# 获取列中的所有字符，并去重
unique_values <- unique(data$column_name)
# 使用factor函数将字符转换为不同的数字
data$column_name <- as.numeric(factor(data$column_name, levels = unique_values))

#########提取

# 从第一列提取基因名称
genename=as.character(DEG[, 1])
genename=as.character(DEG["TRAPPC1", 1])
gene_map=bitr(genename,fromType="SYMBOL",toType="ENTREZID",OrgDb=org.Hs.eg.db)

# # # # # 合并

# 通过符号将gene_map与DEG合并
#即只保留两个数据框中 "SYMBOL" 列具有相同值的行，并将这些行合并成一个新的数据框。
gene_map <- inner_join(gene_map, DEG, by = "SYMBOL")
# 在 clinical_data 中匹配病人ID，获取对应的临床信息（以y为主，可以理解为保留全部，不匹配的会自动填充，若x不够，会按照x的数量走。）
process_TCGA_clinical <- merge(sample_names,  selected_columns, by.x = "patient", by.y = "Patient", all.x = TRUE)
# 组合，要求数据框完全匹配，一般用不到
expr_count = cbind(gene_type=data1$gene_type,gene_name=data1$gene_name,counts)

# # # # # # #去重

#去重方式为取组内最大值，但不一定在同一行，慎用
counts_m <- aggregate( . ~ gene_name,data=counts_m, max) 
#去重duplicated()函数会保留第一个出现的重复行，而将后续的重复行都删除掉。
exp <- exp[!duplicated(exp$id),]
#去重平均值方式去重
counts_m<- ddply(counts_m, .(gene_name), function(x) colMeans(x[, 3:ncol(x)]))
counts_m <- aggregate( . ~ gene_name,data=counts_m, mean)

# 获取排序后的索引
sorted_index <- order(gene_map[, 3], decreasing = TRUE)
# 根据排序后的索引重新排列数据框
gene_map <- gene_map[sorted_index, ]
# 假设您的数据框名为 df，FC 列是需要排序的列
sorted_df <- df[order(abs(df$FC)), ]
#降序排序，可以使用 rev() 函数来反转排序的索引：
sorted_df <- df[rev(order(abs(df$FC))), ]






#8.处理复杂数据
#将聚类结果的基因分列
Go1 <- separate_rows(Go1, genes, sep = "/")
#赋值
#如果x列的某行值小于0，则将该行y列填充为0
lfc$y[lfc$x < 0] <- 0
将数据从宽格式转换为长格式，然后将行名和原来的列名作为新的列
exp1<-tibble::rownames_to_column(exp1, var = "Gene")
exp2 <-pivot_longer(exp1,-Gene)




# 10. 储存表格(更建议第一种方便调用)
saveRDS(clinical,file = paste0(project,"_clinical.rds"))
write.csv(geo$exp,paste0(gse,"_results.csv"))
# 输出图片
filename <- "箱式小提琴"
index <- 1
while (file.exists(paste0(filename, "-", index, ".pdf"))) {
  index <- index + 1
}
pdf(paste0(filename, "-", index, ".pdf"), onefile = FALSE,width =10, height = 7.5)

DEGS$plots

dev.off()
saveRDS('TCGA-COAD_counts_m',file = "TCGA-COAD_counts_m.rds")


pdf(paste0(project,"/",project, "-","Go_Reactomeresult_弦图.pdf"), width = 10, height = 8)
# 使用 chord_dat() 函数创建 chord 图所需的数据
par(
  pin = c(12, 9),  # 设置绘图设备的尺寸
  cex = 0.6,       # 设置全局文字缩放因子
  mar = c(1, 1, 1, 1) # 设置绘图区域的边距
)
# 使用 GOChord() 函数绘制 chord 图，并设置参数(如果报错行数不匹配（KEGG中），只需要把Go4调少就行了，另外颜色数量也要相应调整)
p3<-GOChord(chord3, space = 0.01, gene.order = 'logFC', gene.space = 0.25, gene.size = 3,lfc.col=c(rgb(1, 0, 0, 0.9), rgb(1, 1, 1, 0)), ribbon.col=brewer.pal(7, "Set3"))
print(p3)#这一步可以防止某些文件错误
dev.off()


............................................................
geo初解
# 将表达数据中的Probe ID转换为符号（symbol）此方法随机删除，不可取，仅参考
geo$exp <- trans_array(geo$exp, ids_570)
# 提取符号（symbol），注意symbol大小写
# 当您需要处理和分析结构化数据，并且需要进行各种数据操作和统计分析时
# 创建一个新的数据框ids，只包含SYMBOL列的数据，并去除重复值
head(ids_570$symbol)
str(ids_570)
ids <- as.data.frame(ids_570$symbol)
# 当您需要存储和处理单一类型的一维数据集合，
# 并且进行数学运算、索引、切片、循环迭代等操作时，向量是一个非常有用的数据结构。
ids_v2 <- ids_570$symbol

#正则替换，*+？{2,6}指数量，[a-z][az]列举字符，^指开头$结尾，\b设置边界，\dD数字\wW字母\sS其他键盘符,R中似乎用双杠
#.换行符外的任意字符，   |或，与小括号搭配，<.+> <.+?>判断括号时的长匹配与短匹配，少用
#（）用来分隔为不同部分并从左到右编号，"\\1_\\1_"将编号重组并添加需要的连接符，
gsub("([ab])", "\\1_\\1_", "abc and ABC")
#\\(用来代表真括号，.*常连用，下边用来提取括号里的内容(gsub替换只能替换成一个值不能是数列，也就是说前二个位置只能是字符，当第一个位置代表全名时gsub作用相当于提取)
#gsub：基本功能：替换，附加功能：提取，高级功能：前置判断（可用于避免对同一文件重复进行命名）
a <- gsub(".*\\((.*)\\).tif", "\\1", "001-001-blue(300ms).tif")
#编辑分组文字
group <- ifelse(str_detect(geo$pd$title, "ost"), "ost", "con")
# 三种文字编辑方式根据需要选择，其实不如上面的判断句来的好用
group_list$group <- stringr::str_remove(group_list$group, "hMSC")
group_list$group <- str_replace(group_list$group, ".*-(.*?)_.*", "\\1")
group_list$group <- str_extract(group_list$group, ".*-(.*?)_")
# 如果不是数据框，则尝试将其转换为数据框类型
exp <- as.data.frame(exp)
#如果执行 pdf() 函数后输出的 PDF 文件是空的，可能是因为在调用 pdf() 函数之后没有进行任何绘图操作，或者绘图操作没有被保存到 PDF 文件中。
pdf(paste0("gse","-volcano.pdf"),onefile=FALSE)#输出文件
DEGS$plots
dev.off()
.............................................
#对别人的代码进行更改，记住一定要检查的失误：
参数没有全替换，参数替换了却没有运行赋值，参数大小写失误
。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。
# TCGA

# 加载所需的包
library(data.table)
library(tinyarray)
#数据清洗：data_cleaning.R文件
source("data_cleaning.R")
cleaned_data <- clean_data(counts)
# 从文件中读取LUAD数据并转换为数据框
LUAD <- fread("TCGA-LUAD.htseq_counts.tsv.gz")
LUAD <- as.data.frame(LUAD)
# 将mRNA表达数据的行名称转换为数据框形式
ids_m <- as.data.frame(rownames(counts_m))
ids_m <-ids_m [c(1,1)]
。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。
#WGCNA
source("标量化geo临床信息.R")
processed_pd <- process_clinical_data(geo$pd)
source("通过聚类进行样品筛选.R")
process_sample_data(geo$exp, 260)
#<<- ("在公式中影响环境表格")
processed_geo <<- processed_geo
。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。
#ggplot2
# 控制坐标轴线的外观
theme(axis.line = element_line(color = "black", size = 1, linetype = "solid"))

# 控制刻度线的外观
theme(axis.ticks = element_line(color = "black", size = 0.5, linetype = "dashed"),
      axis.ticks.length = unit(0.2, "cm"))  # 设置刻度线长度

# 控制刻度标签的外观
theme(axis.text = element_text(color = "black", size = 12, angle = 45, hjust = 0.5, vjust = 0.5))

# 控制坐标轴标题的外观
theme(axis.title = element_text(color = "black", size = 14, face = "bold"))

# 添加或移除绘图区域的边框
theme(panel.border = element_rect(color = "black", fill = NA, size = 1))
。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。。
#看富集结果
View(KEGG_gseresult@result)