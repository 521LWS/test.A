
prepare<- function() {
  clinical_data=""
}  





#clinical_data=geo_pd
process_clinical_data <- function(clinical_data) {
  # 根据关键词提取必要的列，并标量化
  keyword <- c("Patient_ID", "Age", "Sex", "gender", "Disease_Status", "diagno")
  selected_columns <- clinical_data[, grepl(paste(keyword, collapse = "|"), names(clinical_data), ignore.case = TRUE)]
  
  selected_columns <- selected_columns[c(1, 2, 3)]
  colnames(selected_columns) <- c("Age", "Disease_Status", "gender")
  
  # 将年龄转换为数值标量使用正则表达式删除非数字字符，并将列转换为数值型数据
  selected_columns$Age <- as.numeric(gsub("[^0-9.]", "", selected_columns$Age))
  # 如果你希望将无效的数值（例如空白字符或无效的数字）替换为特定的值（例如平均值）
  # mean_age <- mean(selected_columns$Age, na.rm = TRUE)
  # df$Age[is.na(selected_columns$Age) | selected_columns$Age == ""] <- mean_age
  # 将疾病状态转换为数值标量（包含Healthy或者是无效值）
  # selected_columns$Disease_Status <- ifelse(grepl("Healthy", selected_columns$Disease_Status, ignore.case = TRUE) | is.na(selected_columns$Disease_Status), 0, 1)
  # 使用factor函数将字符转换为不同的数字(空格也被命名了，如果不删除无效值就不会被命名)
  unique_values <- unique(selected_columns[["Disease_Status"]])
  selected_columns[["Disease_Status"]] <- as.numeric(factor(selected_columns[["Disease_Status"]], levels = unique_values))

  # 将性别转换为二进制编码
  selected_columns$gender <- ifelse(grepl("Female", selected_columns$gender, ignore.case = TRUE), 0, 1)
  
  # 替换 NA
  selected_columns[["Disease_Status"]] <- ifelse(is.na(selected_columns[["Disease_Status"]]), "-1", selected_columns[["Disease_Status"]])
  
  # 现在将数值向量转换为数值类型
  selected_columns[["Disease_Status"]] <- as.numeric(selected_columns[["Disease_Status"]])
  
  
  
  
  return(selected_columns)
}
#调用函数
#source("函数\\函数_标量化geo临床信息.R")
#processed_pd <- process_clinical_data(geo$pd)






















