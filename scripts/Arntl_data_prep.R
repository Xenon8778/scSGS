library(Seurat)
library(dplyr)
library(ggplot2)
source("https://raw.githubusercontent.com/Xenon8778/scSGS/main/R/scType_anno.R")

# Loading Data
tab = readMM('Circadian/data/ZT0/matrix.mtx.gz')
colnames(tab) = read.table('Circadian/data/ZT0/features.tsv.gz')[[2]]
rownames(tab) = read.table('Circadian/data/ZT0/barcodes.tsv.gz')[[1]]
tab = t(tab)
ZT0 = CreateSeuratObject(CreateAssayObject(tab, min.cell = 200, min.features = 3))
ZT0$Batch = 'ZT0'

tab = readMM('Circadian/data/ZT6/matrix.mtx.gz')
colnames(tab) = read.table('Circadian/data/ZT6/features.tsv.gz')[[2]]
rownames(tab) = read.table('Circadian/data/ZT6/barcodes.tsv.gz')[[1]]
tab = t(tab)
ZT6 = CreateSeuratObject(CreateAssayObject(tab, min.cell = 200, min.features = 3))
ZT6$Batch = 'ZT6'

tab = readMM('Circadian/data/ZT12/matrix.mtx.gz')
colnames(tab) = read.table('Circadian/data/ZT12/features.tsv.gz')[[2]]
rownames(tab) = read.table('Circadian/data/ZT12/barcodes.tsv.gz')[[1]]
tab = t(tab)
ZT12 = CreateSeuratObject(CreateAssayObject(tab, min.cell = 200, min.features = 3))
ZT12$Batch = 'ZT12'

tab = readMM('Circadian/data/ZT18/matrix.mtx.gz')
colnames(tab) = read.table('Circadian/data/ZT18/features.tsv.gz')[[2]]
rownames(tab) = read.table('Circadian/data/ZT18/barcodes.tsv.gz')[[1]]
tab = t(tab)
ZT18 = CreateSeuratObject(CreateAssayObject(tab, min.cell = 200, min.features = 3))
ZT18$Batch = 'ZT18'

merged_data = merge(ZT0,c(ZT6,ZT12,ZT18))
so = merged_data

# Filtering
selected_f <- rownames(so)[Matrix::rowSums(so) > 15]
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^mt-")
so.filt <- subset(so,
                  subset = nFeature_RNA > 500 & percent.mt < 5, features = selected_f)
so.filt

# perform visualization and clustering steps
so.filt <- NormalizeData(object = so.filt)
so.filt <- FindVariableFeatures(object = so.filt)
so.filt <- ScaleData(object = so.filt)
so.filt <- RunPCA(object = so.filt)
ElbowPlot(so.filt,ndims = 50)
so.filt <- FindNeighbors(object = so.filt, dims = 1:40)
so.filt <- FindClusters(object = so.filt)
so.filt <- RunUMAP(object = so.filt, dims = 1:40)
p1 = DimPlot(object = so.filt, reduction = "umap", group.by = 'Batch')
p1

# saving
saveRDS(so.filt,'Circadian/data/SMC_annot.rds')

res = FindMarkers(so.annot, ident.1 = 'KO', ident.2 = 'WT', group.by = 'Batch',
                  logfc.threshold = 0.1, min.pct = 0.1) %>%
  mutate(genes = rownames(.))
FeaturePlot(so.filt, features = 'Clock', order = T, split.by = 'Batch')
VlnPlot(so.filt, features = c('Clock','Arntl'), group.by = 'Batch', stack = T,
        flip = T) +
  geom_boxplot(width = 0.2)
table(so.annot$celltype,so.annot$Batch)
