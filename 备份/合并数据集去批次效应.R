

library(limma)
library(stringr)

#1.预处理 
rt1<-GSE31628OE$exp[ ,c(1,2,5,6)]
rt2<-GSE31628OE$exp[ ,c(3,4,7,8)]

rt1<-normalizeBetweenArrays(rt1)
rt2<-normalizeBetweenArrays(rt2)





#2. 合并数据集
#获取共同基因 
gene<-intersect(rownames(rt1), rownames(rt2)) 
#合并
allmerge=cbind(rt1[gene,],rt2[gene,])

allmerge=geo_exp


#3. 根据分组信息去批次
group_list=c(rep("control", 8), rep("si", 8)) 
#批次信息 
batchType=c(rep("GsE1", 1), rep("GsE2", 1),rep("GsE3", 1), rep("GsE4", 1),rep("GsE5", 1), rep("GsE6", 1),rep("GsE7", 1), rep("GsE8", 1),rep("GsE1", 1), rep("GsE2", 1),rep("GsE3", 1), rep("GsE4", 1),rep("GsE5", 1), rep("GsE6", 1),rep("GsE7", 1), rep("GsE8", 1))
#批次信息 
design=model.matrix(~group_list) 
#对数据进行批次矫正，输出矫正后的表达数据 
merge_limma1<-removeBatchEffect(allmerge, batch=batchType, design=design) 
##看一下去除之后的箱线图 
boxplot(merge_limma,outline=FALSE,notch=T,las=2,
        col=c(rep("red",10),rep("blue",10)),
        names=batchType,
        cex.axis=0.7)#出箱线图
geo_exp<-merge_limma[ ,c(1,2,5,6,3,4,7,8)]
boxplot(allmerge,outline=FALSE,notch=T,las=2,
                       col=c(rep("red",10),rep("blue",10)),
                       names=batchType,
                      cex.axis=0.7)#出箱线图
geo_exp=merge_limma1
