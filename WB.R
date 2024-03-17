setwd("D:\\Temp")

# 获取所有以 .scn 结尾的文件
all.file <- list.files(pattern="\\.scn", include.dirs=TRUE, recursive=TRUE, full.names=FALSE)

# 保存旧文件名
old.name <- all.file

# 替换 "总" 为 "t-"
new.name <- gsub(pattern="总", replacement="t-", old.name, ignore.case=TRUE)

# 替换 "Exposure" 为 "Exp"
new.name <- gsub(pattern="Exposure", replacement="Exp", new.name, ignore.case=FALSE)

# 删除日期中的 "-"
new.name <- gsub(pattern="(20[12][0-9])(-)([0-9]{2})(-)([0-9]{2})", replacement="\\1\\3\\5", new.name)

# 连接日期和时间，用 "_"
new.name <- gsub("(?<=20[12][0-9][0-1][0-9][0-9]{2}) (?=[0-9]{2} 时)", "_", new.name, perl=TRUE)

# 将 "时" 替换为 "h"
new.name <- gsub("(?<=[0-9]{2}) 时 (?=[0-9]{2} 分)", "h", new.name, perl=TRUE)

# 将 "分" 替换为 "m"，根据情况添加后缀
new.name<-gsub("(?<=[0-9]{2}) 分(?=_)", "m", new.name, perl=TRUE); #自动曝光
new.name <- gsub("(?<=[0-9]{2}) 分(?= *[0-9])", "m_Man_", new.name, perl=TRUE) # 手动曝光

new.name <- gsub("(?<=[0-9]{2}) 分 *[wW] *(?=.scn)", "m_w", new.name, perl=TRUE) # 白光
new.name <- gsub("(?<=[0-9]{2}) 分 *(?=.scn)", "m", new.name, perl=TRUE) # 白光但是漏了w，我们也不臆断加w

# 删除 "administrator" 字段
new.name <- gsub("administrator ", "", new.name, ignore.case=TRUE, perl=TRUE)

# 添加 () 到文件名
new.name <- gsub("(.*/)(20[123][0-9][0-9]{4}[_0-9a-zA-Z\\. \\+\\-]+)(.tif|.scn)","\\1\\(\\2\\)\\3", new.name, perl=TRUE)

# 将文件夹名添加到文件名中
new.name2 <- gsub("(.*)(/)(.*)(/)(\\(20[123][0-9][0-9]{4}[\\_0-9a-zA-Z\\.\\- \\+]+\\))(.tif|.scn)","\\1\\2\\3\\4\\1 \\3 \\5\\6", new.name, perl=TRUE)

# 连接日期和 gel，分开 gel 和 antibody
new.name3 <- gsub("(.*/)([^ /]+)( )([^ \\./]+)(\\.)([^ \\./]+)( )([^ /]+)(.tif|.scn)","\\1\\2_\\4 \\6\\7\\8\\9", new.name2, perl=TRUE)

# 重命名文件
file.rename(all.file, new.name3)

