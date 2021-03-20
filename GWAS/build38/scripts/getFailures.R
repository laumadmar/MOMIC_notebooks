args <- commandArgs(trailingOnly = TRUE)
data=read.table(args[1],h=F,skip=1)
data$fam=sapply(strsplit(as.character(data$V1), ":"),'[',1)
data$pat=sapply(strsplit(as.character(data$V1), ":"),'[',2)
#data
thres = read.table(args[2],header = TRUE)
failures=subset(data,V2>thres[1,1] | V2<thres[1,2] | V3>thres[1,3] | V3<thres[1,4], select=c(fam,pat))
#failures
write.table(failures,args[3],sep=" ",row.names=FALSE, col.names=FALSE, quote=FALSE)