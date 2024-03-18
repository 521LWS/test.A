

#1.安装两个工具（已安装则忽略）
install.packages("devtools")
install.packages("roxygen2")

#2.在指定路径创建一个空包，这里creat的名字就是包的最终名字
setwd("C:/Users/21989/Desktop/new")
library(devtools)
create("immune.rename")
setwd("C:/Users/21989/Desktop/new/immune.rename")
#3. 
# 在R文件夹中创建一个R文件，输入R语言以及注释这里建议用gpt改，需要包裹成一个函数{}，
#**在R包中，如果函数没有被导出，无法看到和调取（通过在注释的结尾处添加@export）
file.create("R/immune.rename.R")

# 使用roxygen2生成文档到man文件夹：
library(roxygen2)
roxygenize()
#4. 构建和安装包：

devtools::build()
devtools::install()
#5. 测试：
library(immune.rename)
immune.rename("C:/Users/21989/Desktop/new")

