G4=readRDS("/hwfssz1/ST_SUPERCELLS/P22Z10200N0661/haoshijie/FIGURE3/01.analysis/C03638G4/C03638G4.rds")
G4$mt.percent=PercentageFeatureSet(G4,pattern = "^MT-")
G4=subset(G4,subset=mt.percent<15)
G4=subset(G4,subset=scDblFinder.class=="singlet"&nuclei<2)
#G4=subset(G4,subset=nFeature_RNA> (mean(G4$nFeature_RNA)-2*sd(G4$nFeature_RNA)) & nFeature_RNA<(mean(G4$nFeature_RNA)+2*sd(G4$nFeature_RNA)))
#G4=subset(G4,subset=nCount_RNA> (mean(G4$nCount_RNA)-2*sd(G4$nCount_RNA)) & nFeature_RNA<(mean(G4$nCount_RNA)+2*sd(G4$nCount_RNA)))
G4=NormalizeData(G4,verbose = F)
G4=FindVariableFeatures(G4,verbose = F)
G4=ScaleData(G4,verbose = F)
G4=RunPCA(G4,verbose = F)
G4=RunUMAP(G4,dims=1:15,verbose = F)
G4=FindNeighbors(G4,verbose = F,dims=1:15)
G4=FindClusters(G4,verbose = F)

E2=readRDS("/hwfssz1/ST_SUPERCELLS/P22Z10200N0661/haoshijie/FIGURE3/01.analysis/C03640E2/C03640E2.rds")
E2$mt.percent=PercentageFeatureSet(E2,pattern = "^MT-")
E2=subset(E2,subset=mt.percent<15)
E2=subset(E2,subset=scDblFinder.class=="singlet"&nuclei<2)
#E2=subset(E2,subset=nFeature_RNA> (mean(E2$nFeature_RNA)-2*sd(E2$nFeature_RNA)) & nFeature_RNA<(mean(E2$nFeature_RNA)+2*sd(E2$nFeature_RNA)))
#E2=subset(E2,subset=nCount_RNA> (mean(E2$nCount_RNA)-2*sd(E2$nCount_RNA)) & nFeature_RNA<(mean(E2$nCount_RNA)+2*sd(E2$nCount_RNA)))
E2=NormalizeData(E2,verbose = F)
E2=FindVariableFeatures(E2,verbose = F)
E2=ScaleData(E2,verbose = F)
E2=RunPCA(E2,verbose = F)
E2=RunUMAP(E2,dims=1:15,verbose = F)
E2=FindNeighbors(E2,verbose = F,dims=1:15)
E2=FindClusters(E2,verbose = F)

E2$sample="C03640E2"
G4$sample="C03638G4"

C4=readRDS("/hwfssz1/ST_SUPERCELLS/P22Z10200N0661/haoshijie/FIGURE2/01.S1_staining_analysis/B03404C4.rds")
C4$mt.percent=PercentageFeatureSet(C4,pattern = "^MT-")
C4=subset(C4,subset=mt.percent<15)
C4=subset(C4,subset=scDblFinder.class=="singlet"&nuclei<2)
C4=NormalizeData(C4,verbose = F)
C4=FindVariableFeatures(C4,verbose = F)
C4=ScaleData(C4,verbose = F)
C4=RunPCA(C4,verbose = F)
C4=RunUMAP(C4,dims=1:15,verbose = F)
C4=FindNeighbors(C4,verbose = F,dims=1:15)
C4=FindClusters(C4,verbose = F)
C4$sample="B03404C4"

C6=readRDS("/hwfssz1/ST_SUPERCELLS/P22Z10200N0661/haoshijie/FIGURE2/01.S1_staining_analysis/B03404C6.rds")
C6$mt.percent=PercentageFeatureSet(C6,pattern = "^MT-")
C6=subset(C6,subset=mt.percent<15)
C6=subset(C6,subset=scDblFinder.class=="singlet"&nuclei<2)
C6=NormalizeData(C6,verbose = F)
C6=FindVariableFeatures(C6,verbose = F)
C6=ScaleData(C6,verbose = F)
C6=RunPCA(C6,verbose = F)
C6=RunUMAP(C6,dims=1:15,verbose = F)
C6=FindNeighbors(C6,verbose = F,dims=1:15)
C6=FindClusters(C6,verbose = F)
C6$sample="B03404C6"

TenX=readRDS("/hwfssz1/ST_SUPERCELLS/P22Z10200N0661/haoshijie/FIGURE2/01.S1_staining_analysis/human_pbmc_10x_filtered.rds")
TenX$mt.percent=PercentageFeatureSet(TenX,pattern = "^MT-")
TenX=subset(TenX,subset=mt.percent<15)
TenX=NormalizeData(TenX,verbose = F)
TenX=FindVariableFeatures(TenX,verbose = F)
TenX$sample="10X"

obj.list=list(C4,C6,E2,G4,TenX)
intefeatures=SelectIntegrationFeatures(object.list = obj.list,nfeatures = 3000)
options(future.globals.maxSize=100000000000000000)
plan("multicore", workers = 20)
anchors=FindIntegrationAnchors(object.list = obj.list,anchor.features = intefeatures)
plan("multicore", workers = 1)
obj.integrated=IntegrateData(anchors)

obj.integrated@active.assay="integrated"
obj.integrated=ScaleData(obj.integrated,verbose = F)
obj.integrated=RunPCA(obj.integrated,verbose = F)
obj.integrated=RunUMAP(obj.integrated,dims=1:15,verbose = F)
obj.integrated=FindNeighbors(obj.integrated,verbose = F)
obj.integrated=FindClusters(obj.integrated,verbose = F)


############### Remove low quality cells #####################################
obj.integrated=subset(obj.integrated,subset=(! obj.integrated$seurat_clusters %in% c(12,18))|sample=="10X")
obj.integrated@active.assay="integrated"
obj.integrated=ScaleData(obj.integrated,verbose = F)
obj.integrated=RunPCA(obj.integrated,verbose = F)
obj.integrated=RunUMAP(obj.integrated,dims=1:10,verbose = F)
obj.integrated=FindNeighbors(obj.integrated,verbose = F,dims = 1:10)
obj.integrated=FindClusters(obj.integrated,verbose = F)

obj.integrated=subset(obj.integrated,subset=seurat_clusters!=11|sample=="10X")
obj.integrated@active.assay="integrated"
obj.integrated=ScaleData(obj.integrated,verbose = F)
obj.integrated=RunPCA(obj.integrated,verbose = F)
obj.integrated=RunUMAP(obj.integrated,dims=1:10,verbose = F)
obj.integrated=FindNeighbors(obj.integrated,dims=1:10,verbose = F)
obj.integrated=FindClusters(obj.integrated,verbose = F)

obj.integrated$platform="Stereo-Cell"
obj.integrated$platform[obj.integrated$sample=="10X"]="10X"

saveRDS(obj.integrated,file="/hwfssz1/ST_SUPERCELLS/P22Z10200N0661/haoshijie/FIGURE3/01.analysis/mIF.integration/integrated.S1.staining.mIF.rds")
