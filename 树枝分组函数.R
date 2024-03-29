split_and_fill <- function(df = df, chars = chars) {
  
  # Extend chars vector by repeating its last element
  # 将字符向量扩展，重复其最后一个元素
  chars <- c(chars, tail(chars, n = 1))
  
  # Find the row index of the first occurrence of each character
  # 找到每个字符第一次出现的行索引
  first_occurrences <- sapply(chars, function(char) which(df[, 1] == char)[1])
  
  # Split the data frame into parts based on the first occurrences of characters
  # 根据字符的第一次出现将数据框拆分成几部分
  parts <- lapply(seq_along(first_occurrences), function(i) {
    if (i == 1) {
      start_index <- 1
    } else { 
      if (i == length(first_occurrences)) {
        start_index <- first_occurrences[i]-1
      } else {
        start_index <- first_occurrences[i - 1]-1
      }
    }
    
    if (i == length(first_occurrences)) {
      end_index <- nrow(df)
    } else {
      end_index <- first_occurrences[i] - 2
    }
    
    # Extract the part between start and end index
    # 提取开始和结束索引之间的部分
    if (start_index <= end_index) {
      df[start_index:end_index, ]
    } else {
      # If start_index is greater than end_index, return NULL
      # 如果start_index大于end_index，则返回NULL
    }
  })
  
  # Filter out NULL parts
  # 过滤掉NULL部分
  filled_parts <- Filter(Negate(is.null), parts)
  
  # Assign group labels to each part in the third column
  # 为每个部分的第三列分配组标签
  for (i in seq_along(parts)) {
    if (!is.null(parts[[i]])) {
      parts[[i]][, 3] <- paste0("group", i)
    }
  }
  
  # Combine all parts into a single data frame
  # 将所有部分合并成一个数据框
  do.call(rbind, parts)
  
  # Return the list of parts
  # 返回部分列表
  return(parts)
}

# Call the split_and_fill function with an example data frame and characters
# 使用示例数据框和字符调用split_and_fill函数
parts <- split_and_fill(df = as.data.frame(as.phylo(sampleTree2)$edge), chars = c("52", "65", "53"))
df = sampleTree2$edge

