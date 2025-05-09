library(data.table)
mask=fread("B03404C4.cell_mask.csv")
readcount=fread("B03404C4_1725409058149879815/02.count/B03404C4_raw_barcode_gene_exp.txt")
mergecount=merge(readcount,mask,by=c("x","y"))
mergecount=mergecount[,c("x","y","geneIndex","MIDIndex","readCount","cell")]
library(parallel)
cl <- makeCluster(30)
x=mean(aggregate(mergecount$readCount,by=list(mergecount$cell),sum)$x)
y=mean(aggregate(mergecount$geneIndex,by=list(mergecount$cell),FUN=function(xx){length(unique(xx))})$x)
tab=c(x,y)
for(it in seq(1,12,1)){
i=0.8^it
clusterExport(cl=cl,varlist=c("i"))
res=parApply(cl,X=as.matrix(mergecount),FUN=function(xx){sum(sample(0:1,size=xx[5],replace=T,prob=c(1-i,i)))},MARGIN=1)
mergecount$test=res
mergecount=mergecount[mergecount$test!=0,]
x=mean(aggregate(mergecount$test,by=list(mergecount$cell),sum)$x)
y=mean(aggregate(mergecount$geneIndex,by=list(mergecount$cell),FUN=function(xx){length(unique(xx))})$x)
tab=rbind(tab,c(x,y))
}
tab=as.data.frame(tab)
colnames(tab)=c("readCount","geneCount")
saveRDS(tab,"saturation1.rds")
