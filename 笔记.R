
#1. 中文模式
Sys.setlocale(category = "LC_ALL", locale = "en_US.UTF-8")
install.packages("Rcpp")
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

.....................................................................
#变量
a=1:10
dim(a)=c(2,5)




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






