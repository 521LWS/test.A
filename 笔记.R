
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

remove.packages("Rcpp")
install.packages("languageserver", clean = TRUE)
#重新安装
install.packages("languageserver", repos=NULL, type="source")

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
help.search ("BiocManager")

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


install.packages("ctv")


#这段代码是在R中用于检查是否已安装BiocManager包的语句。
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("WGCNA", dependencies = TRUE)
BiocManager::install("GEOquery")
BiocManager::install("limma")
BiocManager::install("affy")
BiocManager::install("tinyarray",dependencies = TRUE)
install.packages("WGCNA",dependencies = TRUE)
install.packages("tinyarray",dependencies = TRUE)
BiocManager::install("clusterProfiler")
install.packages("ComplexHeatmap")
install.packages("org.Rn.eg.db")

install.packages("AnnoProbe",dependencies = TRUE)


library("WGCNA")
library(tibble)
library("tinyarray")
library(stringr)
library(data.table)
library(AnnoProbe)
library(GEOquery)

  ost1 <geo_download("GSE35958")
  ost1_exp <ostlSexp
  ost1 pd <-ost1spd
  ost1_exp <log2(ost1_exp+1)
  find_anno("GPL570",insta11=T)P首选id转旋包
  library(hgu133plus2.db)
  ids_570 <-toTable(hgu133plus2SYMBOL)
  ost1_exp <trans_array(ost1_exp,ids_570)
  group <-ifelse(str_detect(ost1_pdstitle,"ost"),"ost","con")
  broup <as.factor(group)
  ids <-as.data.frame(ids_570Ssymbol)
  idsSv2 <-ids_570Ssymbol
  DEGS <get_deg_all(
    ost1_exp,
    group,
    ids.
    logFC_cutoff 2,
    pvalue_cutoff =0.05,
    adjust =T,
    entriz F,
    scale_before FALSE,
    n_cutoff 3,
    cluster_cols TRUE,
    annotation_legend FALSE,
    show_rownames FALSE,
    legend FALSE,
    lab NA,
    pkg 4,
    symmetry FALSE,
    heat_union TRUE,
    heat_id =1,
    gene_number =200,
    color_volcano c("#2874C5","grey","#f87669")
    PEGSSplots
   