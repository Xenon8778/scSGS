library(Seurat)
library(dplyr)
library(ggplot2)
source("G:/My Drive/Cai lab/Codes/scType_anno.R")

########################
# PBMC 10K
########################
PBMC.10K = Read10X_h5('PBMC/10k_PBMC_3p_nextgem_Chromium_X_filtered_feature_bc_matrix.h5')
PBMC.10K = CreateSeuratObject(CreateAssayObject(PBMC.10K, min.cells = 15,
                                               min.features = 500))

selected_f <- rownames(PBMC.10K)[Matrix::rowSums(PBMC.10K) > 15]
PBMC.10K[["percent.mt"]] <- PercentageFeatureSet(PBMC.10K, pattern = "^mt-")
PBMC.filt <- subset(PBMC.10K,
                    subset = nFeature_RNA > 500 & percent.mt < 5, features = selected_f)

# perform visualization and clustering steps
PBMC.filt <- NormalizeData(object = PBMC.filt)
PBMC.filt <- FindVariableFeatures(object = PBMC.filt)
PBMC.filt <- ScaleData(object = PBMC.filt)
PBMC.filt <- RunPCA(object = PBMC.filt)
ElbowPlot(PBMC.filt,ndims = 50)
PBMC.filt <- FindNeighbors(object = PBMC.filt, dims = 1:40)
PBMC.filt <- FindClusters(object = PBMC.filt)
PBMC.filt <- RunUMAP(object = PBMC.filt, dims = 1:40)
p1 = DimPlot(object = PBMC.filt, reduction = "umap", group.by = 'seurat_clusters')
p1

# Annotating
annotation = Annotate_cells(GetAssayData(PBMC.filt, layer = 'counts'), Scale = T,
                            tissue = 'Immune system', annot_only = T, ndims = 40)
PBMC.annot = PBMC.filt
PBMC.annot$celltype = factor(annotation)
p3 = DimPlot(object = PBMC.annot, reduction = "umap", group.by = 'celltype',
             label = T, label.box = T, repel = T)+NoLegend()
p3
PBMC.annot$celltype = recode_factor(PBMC.annot$celltype,
                                    'Naive CD4+ T cells' = "CD4 T cells",
                                    'Effector CD4+ T cells' = "CD4 T cells",
                                    'Naive CD8+ T cells' = "T cells",
                                    'CD8+ NKT-like cells' = "T cells",
                                    'Naive B cells' = "B cells",
                                    'Memory B cells' = "B cells",
                                    'Pre-B cells' = 'B cells')
p4 = DimPlot(object = PBMC.annot, reduction = "umap", group.by = 'celltype',
             label = T, label.box = T, repel = T)+NoLegend()
p3+p4

# saving
saveRDS(PBMC.annot,'PBMC/PBMC10K_annot.rds')


########################
# PBMC 20K
########################
PBMC.20K = Read10X_h5('PBMC/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5')
PBMC.20K = CreateSeuratObject(CreateAssayObject(PBMC.20K, min.cells = 15,
                                                min.features = 500))

selected_f <- rownames(PBMC.20K)[Matrix::rowSums(PBMC.20K) > 15]
PBMC.20K[["percent.mt"]] <- PercentageFeatureSet(PBMC.20K, pattern = "^mt-")
PBMC.filt <- subset(PBMC.20K,
                    subset = nFeature_RNA > 500 & percent.mt < 5, features = selected_f)

# perform visualization and clustering steps
PBMC.filt <- NormalizeData(object = PBMC.filt)
PBMC.filt <- FindVariableFeatures(object = PBMC.filt)
PBMC.filt <- ScaleData(object = PBMC.filt)
PBMC.filt <- RunPCA(object = PBMC.filt)
ElbowPlot(PBMC.filt,ndims = 50)
PBMC.filt <- FindNeighbors(object = PBMC.filt, dims = 1:40)
PBMC.filt <- FindClusters(object = PBMC.filt)
PBMC.filt <- RunUMAP(object = PBMC.filt, dims = 1:40)
p1 = DimPlot(object = PBMC.filt, reduction = "umap", group.by = 'seurat_clusters')
p1

# Annotating
annotation = Annotate_cells(GetAssayData(PBMC.filt, layer = 'counts'), Scale = T,
                            tissue = 'Immune system', annot_only = T, ndims = 40)
PBMC.annot = PBMC.filt
PBMC.annot$celltype = factor(annotation)
p3 = DimPlot(object = PBMC.annot, reduction = "umap", group.by = 'celltype',
             label = T, label.box = T, repel = T)+NoLegend()
p3
PBMC.annot$celltype = recode_factor(PBMC.annot$celltype,
                                    'Naive CD4+ T cells' = "CD4 T cells",
                                    'Naive CD8+ T cells' = "T cells",
                                    'Effector CD4+ T cells' = "CD4 T cells",
                                    'CD8+ NKT-like cells' = "T cells",
                                    'Naive B cells' = "B cells",
                                    'Memory B cells' = "B cells",
                                    'Pre-B cells' = 'B cells')
p4 = DimPlot(object = PBMC.annot, reduction = "umap", group.by = 'celltype',
             label = T, label.box = T, repel = T)+NoLegend()
p3+p4

# saving
saveRDS(PBMC.annot,'PBMC/PBMC20K_annot.rds')

########################
# PBMC 5K
########################
PBMC.5K = Read10X_h5('PBMC/SC3pv3_GEX_Human_PBMC_filtered_feature_bc_matrix.h5')
PBMC.5K = CreateSeuratObject(CreateAssayObject(PBMC.5K, min.cells = 15,
                                                min.features = 500))

selected_f <- rownames(PBMC.5K)[Matrix::rowSums(PBMC.5K) > 15]
PBMC.5K[["percent.mt"]] <- PercentageFeatureSet(PBMC.5K, pattern = "^mt-")
PBMC.filt <- subset(PBMC.5K,
                    subset = nFeature_RNA > 500 & percent.mt < 5, features = selected_f)

# perform visualization and clustering steps
PBMC.filt <- NormalizeData(object = PBMC.filt)
PBMC.filt <- FindVariableFeatures(object = PBMC.filt)
PBMC.filt <- ScaleData(object = PBMC.filt)
PBMC.filt <- RunPCA(object = PBMC.filt)
ElbowPlot(PBMC.filt,ndims = 50)
PBMC.filt <- FindNeighbors(object = PBMC.filt, dims = 1:40)
PBMC.filt <- FindClusters(object = PBMC.filt)
PBMC.filt <- RunUMAP(object = PBMC.filt, dims = 1:40)
p1 = DimPlot(object = PBMC.filt, reduction = "umap", group.by = 'seurat_clusters')
p1

# Annotating
annotation = Annotate_cells(GetAssayData(PBMC.filt, layer = 'counts'), Scale = T,
                            tissue = 'Immune system', annot_only = T, ndims = 40)
PBMC.annot = PBMC.filt
PBMC.annot$celltype = factor(annotation)
p3 = DimPlot(object = PBMC.annot, reduction = "umap", group.by = 'celltype',
             label = T, label.box = T, repel = T)+NoLegend()
p3
PBMC.annot$celltype = recode_factor(PBMC.annot$celltype,
                                    'Naive CD4+ T cells' = "CD4 T cells",
                                    'Naive CD8+ T cells' = "T cells",
                                    'Memory CD4+ T cells' = "CD4 T cells",
                                    'CD8+ NKT-like cells' = "T cells",
                                    'Naive B cells' = "B cells",
                                    'Pre-B cells' = 'B cells')
p4 = DimPlot(object = PBMC.annot, reduction = "umap", group.by = 'celltype',
             label = T, label.box = T, repel = T)+NoLegend()
p3+p4

# saving
saveRDS(PBMC.annot,'PBMC/PBMC5K_annot.rds')
