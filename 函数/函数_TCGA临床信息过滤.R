

process_TCGA_clinical <- function(TCGA, clinical_data) {
 
  # 根据关键词提取必要的列，并标量化
  keyword <- c("Patient_ID", "patient_barcode", "anatomic_neoplasm_subdivision", "followup", "age_at", "Sex", "gender", "cancer_status", "Disease_Status", "stage_event_pathologic_stage")
  selected_columns <- clinical_data[, grep(paste(keyword, collapse = "|"), names(clinical_data), ignore.case = TRUE)]
  
  selected_columns <- selected_columns[c(1, 2, 5, 6, 7,ncol(selected_columns))]
  colnames(selected_columns) <- c("Patient", "gender", "Age", "tumer","tumer_type", "Disease_Status")
  # 假设 clinical_data 是您的数据框
  column_names = colnames(TCGA)
  # 创建一个包含列名的数据框，并将其命名为 "Sample_Name"
  sample_names <- data.frame(Sample_Name = column_names)
  # 提取样品名称中的病人ID部分
  sample_names$patient <- substr(sample_names$Sample_Name, 1, 12)
  # 在 clinical_data 中匹配病人ID，获取对应的临床信息
 
  process_TCGA_clinical <- merge(sample_names,  selected_columns, by.x = "patient", by.y = "Patient", all.x = TRUE)
  rownames(process_TCGA_clinical)<- process_TCGA_clinical$Sample_Name
   # 删除原始的ID列
  process_TCGA_clinical <-as.data.frame(process_TCGA_clinical[, c(-1,-2)])
  # 使用 for 循环逐行处理
  for (i in seq_len(nrow(process_TCGA_clinical))) {
    # 检查行名是否满足条件（行名的第 14 和 15 个字符是否大于 10）
    if (as.numeric(substring(rownames(process_TCGA_clinical)[i], 14, 15)) > 10) {
      # 如果条件满足，将 Disease_Status 列的值设为 "normal"
      # 将表格中的某一列转换为字符向量
      process_TCGA_clinical$Disease_Status <- as.character(process_TCGA_clinical$Disease_Status)
      
      # 进行赋值操作
      process_TCGA_clinical[i, "Disease_Status"] <- "normal"
      
    }
  }
  
  return(process_TCGA_clinical)
  }




numeric_TCGA_clinical <- function(process_TCGA_clinical) {
  data=process_TCGA_clinical
  # 循环遍历所有列
  column_names <- names(data)
 
  for (col_name in column_names) {
    # 获取列中的所有字符，并去重
    unique_values <- unique(data[[col_name]])
    
    # 使用factor函数将字符转换为不同的数字(空格也被命名了，如果不删除无效值就不会被命名)
    data[[col_name]] <- as.numeric(factor(data[[col_name]], levels = unique_values))
  }
  
  numeric_TCGA_clinical<-data
  return(numeric_TCGA_clinical)
}


#调用函数
#source("函数\\函数_标量化geo临床信息.R")
#process_TCGA_pd <- process_TCGA_clinical(`TCGA-COAD_counts_m`,`TCGA-COAD_pd`)
#numeric_TCGA_clinical <- numeric_TCGA_clinical(`process_TCGA_pd`)