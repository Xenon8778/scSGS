library(Seurat)
library(dplyr)
library(ggplot2)
source("SGS/scType_anno.R")


# RNA Data Prep ####
glioblastoma.rna <- as.sparse(read.csv(file = "Ccr2/RNA/GSE163120_Mouse.GBM.KO1_2_3.WT1_2_3.filtered.gene.bc.matrix.csv.gz",
                                       sep = ",", header = TRUE, row.names = 1, check.names = T))
## Making Seurat Object
glioblastoma.rna <- CreateSeuratObject(CreateAssayObject(glioblastoma.rna))
glioblastoma.rna

# add annotation
annot = read.csv('Ccr2/RNA/GSE163120_annot.Mouse.GBM.KO1_2_3.WT1_2_3.csv.gz',sep = ",", header = T)
annot = annot %>% select(cluster,sample)
glioblastoma.rna$cluster = annot$cluster
glioblastoma.rna$sample = annot$sample

# filtering
selected_f <- rownames(glioblastoma.rna)[Matrix::rowSums(glioblastoma.rna) > 15]
glioblastoma.rna[["percent.mt"]] <- PercentageFeatureSet(glioblastoma.rna, pattern = "^mt-")
glioblastoma.filt <- subset(glioblastoma.rna,
                            subset = nFeature_RNA > 500 & percent.mt < 5, features = selected_f)


# saving
saveRDS(glioblastoma.filt,'Ccr2/RNA/RNA_glioblastoma.rds')

# perform visualization and clustering steps
glioblastoma.filt <- NormalizeData(object = glioblastoma.filt)
glioblastoma.filt <- FindVariableFeatures(object = glioblastoma.filt)
glioblastoma.filt <- ScaleData(object = glioblastoma.filt)
glioblastoma.filt <- RunPCA(object = glioblastoma.filt)
ElbowPlot(glioblastoma.filt,ndims = 50)
glioblastoma.filt <- FindNeighbors(object = glioblastoma.filt, dims = 1:40)
glioblastoma.filt <- FindClusters(object = glioblastoma.filt)
glioblastoma.filt <- RunUMAP(object = glioblastoma.filt, dims = 1:40)
DimPlot(object = glioblastoma.filt, reduction = "umap", group.by = 'seurat_clusters')

# Annotating
annotation = Annotate_cells(GetAssayData(glioblastoma.filt, layer = 'counts'), Scale = T,
                            tissue = 'Immune system', annot_only = T, ndims = 40)
glioblastoma.annot = glioblastoma.filt
glioblastoma.annot$celltype = factor(annotation)
DimPlot(object = glioblastoma.annot, reduction = "umap", group.by = 'celltype')

# saving
saveRDS(glioblastoma.annot,'Ccr2/RNA/RNA_glioblastoma_UMAP.rds')

# CITE Data Prep ####
count.rna <- as.sparse(read.csv(file = "Ccr2/CITE/GSE163120_Citeseq_Mouse.GBM.KO4_WT4.filtered.RNA.feature.bc.matrix.csv.gz",
                               sep = ",", header = TRUE, row.names = 1))

glioblastoma.rna <- CollapseSpeciesExpressionMatrix(count.rna)

glioblastoma.adt <- as.sparse(read.csv(file = "Ccr2/CITE/GSE163120_Citeseq_Mouse.GBM.KO4_WT4.filtered.ADT.feature.bc.matrix.csv.gz",
                               sep = ",", header = TRUE, row.names = 1))

all.equal(colnames(glioblastoma.rna), colnames(glioblastoma.adt))


## Making Seurat Object
glioblastoma <- CreateSeuratObject(glioblastoma.rna)
### filtering
selected_f <- rownames(glioblastoma)[Matrix::rowSums(glioblastoma) > 15]
glioblastoma <- glioblastoma[selected_f,]
Assays(glioblastoma)

# create a new assay to store ADT information
adt_assay <- CreateAssay5Object(counts = glioblastoma.adt)

# add this assay to the previously created Seurat object
glioblastoma[["ADT"]] <- adt_assay

Assays(glioblastoma)

# add annotation
annot = read.csv('Ccr2/CITE/GSE163120_annot.Citeseq_Mouse.GBM.KO4_WT4.csv.gz',sep = ",", header = F)[,-c(7,8)]
colnames(annot) = annot[1,]
annot = annot[-1,] %>% select(cluster,sample)
glioblastoma$cluster = annot$cluster
glioblastoma$sample = annot$sample

## filtering
glioblastoma.filt <- subset(glioblastoma, subset = nFeature_RNA > 500)
Assays(glioblastoma.filt)
dim(glioblastoma.filt)

# perform visualization and clustering steps
DefaultAssay(glioblastoma.filt) <- "RNA"
glioblastoma.filt <- NormalizeData(glioblastoma.filt)
glioblastoma.filt <- FindVariableFeatures(glioblastoma.filt)
glioblastoma.filt <- ScaleData(glioblastoma.filt)
glioblastoma.filt <- RunPCA(glioblastoma.filt, verbose = FALSE)
ElbowPlot(glioblastoma.filt)
glioblastoma.filt <- FindNeighbors(glioblastoma.filt, dims = 1:40)
glioblastoma.filt <- FindClusters(glioblastoma.filt, resolution = 0.8, verbose = FALSE)
glioblastoma.filt <- RunUMAP(glioblastoma.filt, dims = 1:40)
DimPlot(glioblastoma.filt, label = TRUE, reduction = 'umap', group.by = 'cluster')

# Normalize ADT data
glioblastoma.filt <- NormalizeData(glioblastoma.filt, normalization.method = "CLR", margin = 2, assay = "ADT")

# saving
saveRDS(glioblastoma.filt,'Ccr2/CITE/CITE_glioblastoma.rds')

