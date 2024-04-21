prepare<- function() {
  TCGA=`TCGA-BRCAprocessed_counts_m`
  clinical_data=`TCGA-BRCA_pd`
 
  
  
}  

process_TCGA_clinical <- function(TCGA, clinical_data) {

 
  # 根据关键词提取必要的列，并标量化
  keyword <- c( "patient_barcode", "person_neoplasm_cancer_status", "vital_status","days_to_last_followup","days_to_death","relative_family_cancer_history","history_hepato_carcinoma_risk_factors", "age_at_initial_pathologic_diagnosis","histological_type","stage_event_tnm_categories", "gender","height" ,"weight", "neoplasm_histologic_grade","stage_event_pathologic_stage")
  selected_columns <- clinical_data[, grep(paste(keyword, collapse = "|"), names(clinical_data), ignore.case = TRUE)]
  library(dplyr)
  selected_columns$event <- case_when(
    stringr::str_detect(selected_columns$days_to_death, "\\d") ~ 1,
    stringr::str_detect(selected_columns$days_to_last_followup, "\\d") ~ 0,
    stringr::str_detect(selected_columns$vital_status, "Dead") ~ 1,
    stringr::str_detect(selected_columns$vital_status, "Alive") ~ 0,
    TRUE ~ NA_real_
  )
  
 
  
  selected_columns$time<- case_when(
    stringr::str_detect(selected_columns$days_to_death, "\\d") ~ as.numeric(selected_columns$days_to_death),
    TRUE ~ NA_real_
  )

    # 使用正则表达式从Category列中提取T、N、M的内容
  selected_columns $T <-as.numeric(sub(".*T(\\d).*", "\\1",  selected_columns $stage_event_tnm_categories)) 
  selected_columns $N <-as.numeric(sub(".*N(\\d).*", "\\1", selected_columns $stage_event_tnm_categories)) 
  selected_columns $M <- as.numeric(sub(".*M(\\d).*", "\\1", selected_columns $stage_event_tnm_categories)) 

  if(!is.null(selected_columns$weight)){selected_columns$BMI <- as.numeric(selected_columns$weight) / as.numeric(selected_columns$height)
  }
  

  colname1=c("bcr_patient_barcode","time","event","gender","age_at_initial_pathologic_diagnosis","histological_type","stage_event_pathologic_stage","T","N","M", "person_neoplasm_cancer_status","relative_family_cancer_history","history_hepato_carcinoma_risk_factors","BMI")
  # 包含关键词的列名
  selected_columns1 <-  selected_columns[, grep(paste(colname1, collapse = "|"),names(selected_columns))]
  
  # 按照指定顺序排序
  selected_columns <- selected_columns[ ,order(match(colnames(selected_columns),colname1))]
  colnames(selected_columns)[1:10]<-c("Patient","time","event","gender","Age","histological_type","Disease_Status","T","N","M")
  
  
  
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
  
  process_TCGA_clinical=process_TCGA_clinical[order(process_TCGA_clinical[,"Disease_Status"], decreasing = FALSE),]
  data=process_TCGA_clinical
  #将空格换为无效值
  library (dplyr)
  data <- mutate_all(data, as.character)
  data<- data %>% mutate_all (na_if,"")
  data1<- mutate_all(data, as.numeric)
  
  # 循环遍历每一列
  for (i in 1:ncol(data)) {
    # 获取列中的所有字符，并去重
    unique_values <- unique(data[[i]])
    # 使用factor函数将字符转换为不同的数字
    data[[i]] <- as.numeric(factor(data[[i]], levels = unique_values))
  }
  
  data<- mutate_all(data, as.numeric)
  data$Age<-data1$Age
  data$time<-data1$time
 
  if(is.null(data$BMI)){
    data$event<-data1$event
  } else {
    data$BMI<-data1$BMI
  }
    
  numeric_TCGA_clinical<-data
  return(numeric_TCGA_clinical)
}

#调用函数
#source("函数\\函数_标量化geo临床信息.R")
#process_TCGA_pd <- process_TCGA_clinical(`TCGA-COAD_counts_m`,`TCGA-COAD_pd`)
#numeric_TCGA_clinical <- numeric_TCGA_clinical(`process_TCGA_pd`)