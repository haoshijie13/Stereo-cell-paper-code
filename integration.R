library(Seurat)
obj=readRDS("../Y00309KA/Y00309KA.rds")
obj=NormalizeData(obj)
obj=FindVariableFeatures(obj)
obj=ScaleData(obj)
obj=RunPCA(obj)
obj=FindNeighbors(obj,dims=1:15)
obj=FindClusters(obj)
obj=RunUMAP(obj,dims=1:15)
obj1=readRDS("../Y00309KB/Y00309KB.rds")
obj1=NormalizeData(obj1)
obj1=FindVariableFeatures(obj1)
obj1=ScaleData(obj1)
obj1=RunPCA(obj1)
obj1=RunUMAP(obj1,dims=1:15)
obj1=FindNeighbors(obj1,dims=1:15)
obj1=FindClusters(obj1)

obj$sample='KA'
obj1$sample='KB'
library(future)
plan("multicore", workers = 20)
anchors <- FindIntegrationAnchors(object.list = list(obj,obj1))
plan("multicore", workers = 1)
integrated.obj=IntegrateData(anchorset = anchors)
# Run the standard workflow for visualization and clustering
integrated.obj <- ScaleData(integrated.obj, verbose = FALSE)
integrated.obj <- RunPCA(integrated.obj, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
integrated.obj <- RunUMAP( integrated.obj, reduction = "pca", dims = 1:20)
integrated.obj <- FindNeighbors( integrated.obj, reduction = "pca", dims = 1:20)
integrated.obj <- FindClusters( integrated.obj, resolution = 0.5)

s.gene=apply(as.matrix(cc.genes$s.genes),1,function(xx){gsub("^([[:alpha:]])(.*)", "\\1\\L\\2", xx, perl=TRUE)})
g2m.gene=apply(as.matrix(cc.genes$g2m.genes),1,function(xx){gsub("^([[:alpha:]])(.*)", "\\1\\L\\2", xx, perl=TRUE)})
integrated.obj@active.assay="RNA"
integrated.obj=CellCycleScoring(integrated.obj,s.features=s.gene,g2m.features=g2m.gene)

saveRDS(integrated.obj,"integrated.obj.rds")
