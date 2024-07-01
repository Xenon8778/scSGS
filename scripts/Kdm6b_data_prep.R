library(Seurat)
library(dplyr)
library(ggplot2)
source("SGS/scType_anno.R")

# Data Prep ####
Olig.WT = Read10X('Kdm6b/data/WT/')
Olig.WT = CreateSeuratObject(CreateAssayObject(Olig.WT))
Olig.WT$Batch = 'WT'
Olig.KO = Read10X('Kdm6b/data/KO/')
Olig.KO = CreateSeuratObject(CreateAssayObject(Olig.KO))
Olig.KO$Batch = 'KO'

# Merging WT and KO
merged_data = merge(Olig.WT,Olig.KO,add.cell.ids = c("WT", "KO"))

# filtering
selected_f <- rownames(merged_data)[Matrix::rowSums(merged_data) > 15]
merged_data[["percent.mt"]] <- PercentageFeatureSet(merged_data, pattern = "^mt-")
Olig.filt <- subset(merged_data,
                            subset = nFeature_RNA > 500 & percent.mt < 5, features = selected_f)


# perform visualization and clustering steps
Olig.filt <- NormalizeData(object = Olig.filt)
Olig.filt <- FindVariableFeatures(object = Olig.filt)
Olig.filt <- ScaleData(object = Olig.filt)
Olig.filt <- RunPCA(object = Olig.filt)
ElbowPlot(Olig.filt,ndims = 50)
Olig.filt <- FindNeighbors(object = Olig.filt, dims = 1:40)
Olig.filt <- FindClusters(object = Olig.filt)
Olig.filt <- RunUMAP(object = Olig.filt, dims = 1:40)
p1 = DimPlot(object = Olig.filt, reduction = "umap", group.by = 'seurat_clusters')
p2 = DimPlot(object = Olig.filt, reduction = "umap", group.by = 'Batch')
p1+p2

# Annotating
annotation = Annotate_cells(GetAssayData(Olig.filt, layer = 'counts'), Scale = T,
                            tissue = 'Brain', annot_only = T, ndims = 40,
                            filter_conf = F, res = 1)
Olig.annot = Olig.filt
Olig.annot$celltype = factor(annotation)
#Olig.annot$celltype = recode_factor(Olig.annot$celltype)

p3 = DimPlot(object = Olig.annot, reduction = "umap", group.by = 'celltype')
p2+p3

table(Olig.annot$Batch,Olig.annot$celltype)

# saving
saveRDS(Olig.annot,'Kdm6b/Kdm6b_Neuron.rds')

FeaturePlot(Olig.annot, 'Kdm6b')
