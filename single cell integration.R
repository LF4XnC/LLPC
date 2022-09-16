##step.1 create Seurat objets for each sample
seq.data<-Read10X(path)
obj <- CreateSeuratObject(counts = seq.data, project = project, min.cells = 1)

##step.2 merge all samples 
llpc <- merge(obj1,y=c(obj2,obj3,objn),add.cell.ids=c('sample1','sample2','sample3','samplen'),project=project)

##step.3 remove all ig transcripts
llpc@assays$RNA@counts<-llpc@assays$RNA@counts[rownames(llpc@assays$RNA@counts) %in% c(grep('Ig[hkl][vdjcaemg][1-9]*',rownames(llpc),value = T,invert = T)),]

##step.4 basic single cell analysis  pipline
llpc <- CreateSeuratObject(counts = llpc@assays$RNA@counts, min.cells = 1)
llpc[["percent.mt"]] <- PercentageFeatureSet(llpc, pattern = "^mt-")
VlnPlot(llpc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0)
plot1 <- FeatureScatter(llpc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(llpc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
llpc <- NormalizeData(llpc, normalization.method = "LogNormalize", scale.factor = 10000)
llpc <- FindVariableFeatures(llpc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(llpc), 10)
plot1 <- VariableFeaturePlot(llpc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(llpc)
llpc <- ScaleData(llpc, features = all.genes)
llpc <- RunPCA(llpc, features = VariableFeatures(object = llpc))
print(llpc[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(llpc, dims = 1:2, reduction = "pca")
DimPlot(llpc, reduction = "pca")
DimHeatmap(llpc, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(llpc, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(llpc)
llpc <- FindNeighbors(llpc, dims = 1:15)
llpc <- FindClusters(llpc, resolution = 0.5)
llpc <- RunUMAP(llpc, dims = 1:15, reduction = "pca")
DimPlot(llpc, label = TRUE,reduction = 'umap',group.by = 'orig.ident')
DimPlot(llpc, label = TRUE,reduction = 'umap')
FeaturePlot(llpc, features = c('Bach2','Xbp1','Sdc1','Prdm1','nCount_RNA','percent.mt'),
            ncol = 3, pt.size = 0.5,order = T)

##step.5 exclude cluster with abnormal gene expression 
llpc <- subset(llpc,ident=c(8,10,11),invert=T)#then re-run the pepline above

##step.6 exclude cells with no BCR of heavy chain detected and BCR doublets(more than one pair BCR detected)
llpc$vgene <- ig$v_gene.x
llpc$ctype <- factor(ig$c_type,levels = c('IGHA','IGHG','IGHM','IGHD','IGHE','None'))
llpc$aa33 <- ig$aa33.x
llpc$igl <- ig$locus.y
llpc <- subset(llpc, subset = ctype == 'double',invert=T) %>% subset(subset = ctype == 'ND',invert=T) 
#exlude cells in BCR object
ig <- ig[rownames(ig)%in% colnames(llpc),]
ig <- ig[colnames(llpc),,drop=F]

##step.7-1 for reference dataset, perform integration to correct batch effect
llpc.list <- SplitObject(llpc, split.by = "exp")
llpc.list <- llpc.list[c("e1", "e2", "e3", "e4")]
for (i in 1:length(llpc.list)) {
  llpc.list[[i]] <- NormalizeData(llpc.list[[i]], verbose = FALSE)
  llpc.list[[i]] <- FindVariableFeatures(llpc.list[[i]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
}
llpc.anchors <- FindIntegrationAnchors(object.list = llpc.list, dims = 1:30)
llpc.integrated <- IntegrateData(anchorset = llpc.anchors, dims = 1:30)
#switch to integrated assay. The variable features of this assay are automatically set during IntegrateData
DefaultAssay(llpc.integrated) <- "integrated"
#Run the standard workflow for visualization and clustering
llpc.integrated <- ScaleData(llpc.integrated, verbose = FALSE)
llpc.integrated <- RunPCA(llpc.integrated, npcs = 30, verbose = FALSE)
ElbowPlot(llpc.integrated,ndims = 30)
llpc.integrated <- FindNeighbors(llpc.integrated, dims = 1:20)
llpc.integrated <- FindClusters(llpc.integrated, resolution = 1)
llpc.integrated <- RunUMAP(llpc.integrated, reduction = "pca", dims = 1:20)

##step.7-2 for other dataset, perform integration to predict cluster 
#llpc.query as the dataset, llpc.integrated as the reference dataset after integration
llpc.anchors <- FindTransferAnchors(reference = llpc.integrated, query = llpc.query,
                                        dims = 1:30)
predictions <- TransferData(anchorset = llpc.anchors, refdata = llpc.integrated$cluster, 
                            dims = 1:30)
llpc.query <- AddMetaData(llpc.query, metadata = predictions)

##step 8 futher analysis and visualization








