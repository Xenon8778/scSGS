library(Seurat)
library(dplyr)
library(ggplot2)
library(scSGS)
library(SeuratWrappers)
source("SGS/scType_anno.R")

# RNA Data Prep ####

## Control
control.count <- Read10X('Twocondition/GSE235394_Control')
control <- CreateSeuratObject(CreateAssayObject(control.count, min.cells=3,
                                                min.features=200))
control$sample <- 'Control'

## Cancer
cancer.count <- Read10X('Twocondition/GSE235394_B16F10_metastatic_lung')
cancer <- CreateSeuratObject(CreateAssayObject(cancer.count, min.cells=3,
                                               min.features=200))
cancer$sample <- 'Cancer'

so <- merge(control,cancer)

# filtering
selected_f <- rownames(so)[Matrix::rowSums(so) > 15]
so[["percent.mt"]] <- PercentageFeatureSet(so, pattern="^mt-")
so.filtered <- subset(so, subset=nFeature_RNA > 500 & percent.mt < 5,
             features = selected_f)
so.filtered

# saving
saveRDS(so.filtered, 'Twocondition/Lung_cancer.rds')

# perform visualization and clustering steps
so.filtered <- readRDS('Twocondition/Lung_cancer.rds')
so.filtered <- NormalizeData(object = so.filtered)
so.filtered <- FindVariableFeatures(object = so.filtered)
so.filtered <- ScaleData(object = so.filtered)
so.filtered <- RunPCA(object = so.filtered)
ElbowPlot(so.filtered,ndims = 50)
so.filtered <- FindNeighbors(object = so.filtered, dims = 1:40)
so.filtered <- FindClusters(object = so.filtered, resolution = 0.5)
so.filtered <- RunUMAP(object = so.filtered, dims = 1:40)
DimPlot(object = so.filtered, reduction = "umap", group.by = 'sample')
DimPlot(object = so.filtered, reduction = "umap", group.by = 'seurat_clusters',
        label = T)
FeaturePlot(so.filtered, 'Pim3')

# general capillary (gCap)
FindMarkers(so.filtered, ident.1 = 5, only.pos = T) %>% head(10)
FeaturePlot(so.filtered, c('Plvap','Car4','Amigo2','Col3a1'), alpha = 0.3,
            min.cutoff = 1, ncol = 2)

so.annot <- so.filtered
so.annot$celltype <- so.annot$seurat_clusters
so.annot$celltype <- recode_factor(so.annot$celltype, '0' = 'gCap',
                                      '1' = 'gCap',
                                      '2' = 'gCap',
                                      '3' = 'aCap',
                                      '4' = 'aCap',
                                      '5' = 'gCap',
                                      '6' = 'Vein',
                                      '7' = 'gCap',
                                      '8' = 'gCap',
                                      '9' = 'gCap',
                                      '10' = 'gCap',
                                      '11' = 'gCap',
                                      '12' = 'Unknown',
                                      '13' = 'Fibroblast',
                                      '14' = 'Unknown',
                                      '15' = 'Unknown',
                                      '16' = 'Unknown',
                                      '17' = 'gCap',
                                      '18' = 'Unknown')

DimPlot(object = so.annot, reduction = "umap", group.by = 'celltype')

# saving
saveRDS(so.annot,'Twocondition/Lung_cancer_UMAP.rds')

