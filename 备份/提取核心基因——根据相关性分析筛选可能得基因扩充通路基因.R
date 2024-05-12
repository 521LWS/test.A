#获取核心基因
Go_gseresult_core<-as.data.frame(Go_gseresult)
Go_gseresult_core<-Go_gseresult_core[c("ID","core_enrichment")]
#一个搜索函数
get_second_column <- function(keyword, data_frame) {
  row <- which(data_frame[, 1] == keyword)
  if (length(row) == 0) {
    return(NULL)
  } else {
    core1<-unlist(strsplit(data_frame[row, 2], "/"))
    core1<-as.data.frame(core1)
    return(core1)
  }
}