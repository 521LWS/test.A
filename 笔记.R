
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

#检测Rtool安装与否
system('where make')
#PACKAGE路径
.libPaths()
#程序运行路径
getwd()
setwd()
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
#从github安装
install.packages("devtools")
devtools::install_github("Bioconductor/GO.db")



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
geo初解
#移除环境文件
rm(ids, ids_570, geo_pd, geo_exp, fit, fit2, df, data, merged_data,probe_annotation,genes_expr,geo,exp,exp_exp,ensembl,averaged_data,contrast.matrix,DEG,DEGS)

rm(list = ls())  #删除所有对象
# 将表达数据中的Probe ID转换为符号（symbol）
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
group_list$group <-stringr:str_remove(group_listSgroup,"tissue type:")
# 三种文字编辑方式根据需要选择，其实不如上面的判断句来的好用
group_list$group <- stringr::str_remove(group_list$group, "hMSC")
group_list$group <- str_replace(group_list$group, ".*-(.*?)_.*", "\\1")
group_list$group <- str_extract(group_list$group, ".*-(.*?)_")
# 如果不是数据框，则尝试将其转换为数据框类型
exp <- as.data.frame(exp)
# 将结果转换为 tibble 格式(行名入表命名为gene)
DEG <- tibble::rownames_to_column(DEG, var = "Gene")
# 将第一列转换为行名
rownames(geo_exp)<- geo_exp$symbol
# 删除原始的ID列
geo_exp <- geo_exp[, -1]
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
# 将行名称设置为LUADSEnsemb1_ID列中的值，为了方便其他包处理
row.names(LUAD) <- LUADSEnsemb1_ID
# 去除第一列（假设第一列是行索引）
LUAD <- LUAD[, -1]
