
# 定义函数来处理文件夹中的文件
rename_files <- function(folder_path) {
  # 获取文件列表
  files <- list.files(folder_path, full.names = TRUE)
  
  # 过滤非隐藏文件
  files <- files[!file.info(files)$isdir]
  
  # 遍历文件列表
  for (i in seq_along(files)) {
    file_path <- files[i]
    # 提取文件序号
    n <- i
    # 新文件名
    new_name <- ifelse(n %% 2 == 0, paste0(n / 2, "-green"), paste0((n + 1) / 2, "-blue"))
    # 提取文件类型
    extension <- tools::file_ext(file_path)
    # 构建新文件路径
    new_file_path <- file.path(dirname(file_path), paste0(new_name, ".", extension))
    # 重命名文件
    file.rename(file_path, new_file_path)
  }
}

# 递归函数获取所有子文件夹并运行文件重命名函数
recursive_rename <- function(folder_path) {
  # 获取当前文件夹下的所有子文件夹
  subdirs <- list.dirs(folder_path, recursive = FALSE)
  
  # 如果没有子文件夹，则直接对当前文件夹下的文件进行重命名
  if (length(subdirs) == 0) {
    rename_files(folder_path)
  } else {
    # 否则，对每个子文件夹递归调用该函数
    for (subdir in subdirs) {
      recursive_rename(subdir)
    }
  }
}
# 调用递归函数开始处理
recursive_rename("C:/Users/21989/Desktop/new")
print(new_name)


