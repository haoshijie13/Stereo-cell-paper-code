args=commandArgs(T)
library(data.table)
gem=fread(args[1])
prefix=sub(".gem","",args[1])
dapi=fread(args[2])

#过滤一下DAPI，太小或者太大的DAPI信号是假信号
dapi=dapi[dapi$value %in% names(table(dapi$value)[ table(dapi$value) > 100 & table(dapi$value) <1500]),]

mask=read.csv(args[3],row.names=1)
mask=as.data.table(mask)
measure=read.csv(args[4],row.names=1)
#gem1=gem[gem$label!=0,]
mer=merge(gem,mask,by.x=c("x","y"),by.y=c("x","y"),all=T)
mer=merge(mer,dapi,by.x=c("x","y"),by.y=c("x","y"),all=T)
mer1=mer[!is.na(mer$cell),]
mer1=mer1[!is.na(mer1$geneID),]
nuclei=unique(mer1[,c("cell","value")])
nuclei=nuclei[!is.na(nuclei$value),]
cell=1:length(unique(mer1$cell))
names(cell)=unique(mer1$cell)
gene=1:length(unique(mer1$geneID))
names(gene)=unique(mer1$geneID)
library(Matrix)
print("OK")
mat=sparseMatrix(i=gene[ mer1$geneID],j=cell[ as.character(mer1$cell) ],x=mer1$MIDCount)
rownames(mat)=names(gene)
colnames(mat)=paste0("C",names(cell))
measure$cell=round(measure$cell,digits=0)
rownames(measure)=paste0("C",measure$cell)
library(Seurat)
obj=CreateSeuratObject(mat)
humangene=rownames(obj@assays$RNA@counts)[ grep("^GRCh38",rownames(obj@assays$RNA@counts)) ]
mousegene=rownames(obj@assays$RNA@counts)[ grep("^mm10",rownames(obj@assays$RNA@counts)) ]
obj$human_counts=colSums(obj@assays$RNA@counts[humangene,])
obj$mouse_counts=colSums(obj@assays$RNA@counts[mousegene,])
obj@meta.data$type="undefine"
obj@meta.data$percent=obj$human_counts/(obj$nCount_RNA)
obj@meta.data$type[obj$percent>0.9]="human"
obj@meta.data$type[obj$percent<0.2]="mouse"
obj@meta.data=cbind(obj@meta.data,measure[rownames(obj@meta.data),])

library(scDblFinder)
obj1 <- as.SingleCellExperiment(obj)
obj1 <- scDblFinder(obj1, dbr = 0.1)
obj1 <- as.Seurat(obj1)
obj$scDblFinder.class=obj1$scDblFinder.class

nuclei$cell=paste0("C",nuclei$cell)
obj$nuclei=table(nuclei$cell)[rownames(obj@meta.data)]
obj$nuclei[is.na(obj$nuclei)]=0
obj$nuclei=as.numeric(obj$nuclei)

library(ggplot2)
pdf(paste0(prefix,".pdf"))
ggplot(obj@meta.data,aes(x=mouse_counts,y=human_counts,color=type ))+geom_point()+theme_classic()+scale_color_manual(values=c("red","blue","gray"))+labs(color="")
dev.off()
saveRDS(obj,paste0(prefix,".rds"))
