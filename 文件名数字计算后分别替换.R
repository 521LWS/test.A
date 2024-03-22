setwd("C:/Users/21989/Desktop/new")

# Load the stringr package
library(stringr)

# 获取所有以 .scn 结尾的文件
all.file <- list.files(pattern = "\\.tif", full.names = TRUE, recursive = TRUE)
old.names <- all.file

# 定义一个函数用于替换
replace_function <- function(match) {
  time_in_ms <- as.integer(sub("ms", "", match)) * 16
  sprintf("%ds%03dms", time_in_ms %/% 1000, time_in_ms %% 1000)
}

# 生成新文件名
new.names <- str_replace_all(old.names, "(\\d+)ms", replace_function)

# Print the old and new file names
print(paste("Old Names:", old.names))
print(paste("New Names:", new.names))

# Rename files
file.rename(old.names, new.names)

