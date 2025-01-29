# single-cell Stochastic Gene Silencing (scSGS)
Read the full article [**Here!**](https://doi.org/10.1038/s42003-025-07530-0)
This is the main repository for scSGS - A method to study gene function using the natural stochasticity in gene expression.
The scripts for preprocessing and analyzing publicly available datasets are in the '/scripts' folder.

<img src="https://github.com/Xenon8778/scSGS/assets/61325229/9c629e1e-4b34-456d-a80e-7476622ff6b4" alt="Sample Image" width="75%" height="auto">

# Abstract
Gene expression is a dynamic and stochastic process characterized by transcriptional bursting followed by periods of silence. Single-cell RNA sequencing (scRNA-seq) is a powerful tool to measure transcriptional bursting and silencing at the individual cell level. In this study, we introduce the single-cell Stochastic Gene Silencing (scSGS) method, which leverages the natural variability in single-cell gene expression to decipher gene function. For a target gene g under investigation, scSGS classifies cells into transcriptionally active (g+) and silenced (g-) samples. It then compares these cell samples to identify differentially expressed genes, referred to as SGS-responsive genes, which are used to infer the function of the target gene g. Analysis of real data demonstrates that scSGS can reveal regulatory relationships up- and downstream of target genes, circumventing the survivorship bias that often affects gene knockout and perturbation studies. scSGS thus offers an efficient approach for gene function prediction, with significant potential to reduce the use of genetically modified animals in gene function research.

## Citation
Gupta, S., Cai, J.J. Gene function revealed at the moment of sitochastic gene silencing. Commun Biol 8, 88 (2025). https://doi.org/10.1038/s42003-025-07530-0

# Install scSGS
scSGS requires [presto](https://github.com/immunogenomics/presto) for fast computation of DE genes using the wilcoxon rank-sum test.
```R
if(!require(remotes)) install.packages('remotes')
remotes::install_github("immunogenomics/presto")
remotes::install_github("Xenon8778/scSGS")
```

# Quick Start
### Identification of HVGs
Identify HVGs using the `HVG_splinefit` function. This method leverages a spline-based approach to account for non-linear relationships between mean expression, coefficient of variation, and dropout rate.
```R
HVG_res = HVG_splinefit(GetAssayData(seuratObj, layer='counts'), diptest=TRUE,
                        dropout.filter=TRUE, nHVGs=500)
```                        
### scSGS analysis
scSGS analysis should be performed for highly variable genes only. The input for the algorithm is a  scRNA-seq geneexpression matrix with genes as rows and cells as columns. The output of the algorithm is data.frame find the SGS-responsive genes and their P-values. HVG selection is done using HVG_splinefit which uses each gene's mean expression, CV, and dropout rate to compute variability.
```R
SGS_res <- scSGS(GetAssayData(seuratObj, layer='counts'), GoI = 'nameOfGene', calcHVG = T)
```
# Tutorial
### Import Libraries
Load the essential libraries for this analysis: 'Seurat' for single-cell data analysis and 'scSGS' for the specific functionalities we'll be using.
```{r Loading libraries, message=FALSE}
library(Seurat)
library(scSGS)
```

### Load datasets
For this tutorial, We'll work with a publicly available Peripheral Blood Mononuclear Cells (PBMC) dataset from 10X Genomics.  The raw data can be found [here] (https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz). Read the 10X-formatted data and create a Seurat object. 
```{r Load Seurat data, warning=FALSE}
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "../../scSGS_Extra/data/filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(CreateAssayObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200))
pbmc
```

### Filtering 
To improve data quality, filter out cells with extremely low or high feature counts, which often indicate low-quality or technical noise. We filter cells that have unique feature counts over 2,500 or less than 200.
```{r Filtering}
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
pbmc
```

### Identification of HVGs
Identify HVGs using the `HVG_splinefit` function. This method leverages a spline-based approach to account for non-linear relationships between mean expression, coefficient of variation, and dropout rate. scSGS needs a scRNA-seq Matrix or dgCMatrix object with genes as rows and cells as columns, hence Assay Data needs to be extracted to be used as input.
```{r}
HVG_res = HVG_splinefit(GetAssayData(pbmc,layer = 'counts'), diptest=TRUE,
                        dropout.filter=TRUE, nHVGs=500)

# HVG results
head(HVG_res[HVG_res$HVG == TRUE,])
```
```{r TFs}
# Highly Variable Genes that are known TFs 
HVG_res[(HVG_res$HVG == TRUE & HVG_res$TF == TRUE),][1:5,]
```
### scSGS analysis
Perform scSGS analysis on the selected HVGs using the `scSGS` function. This function identifies genes that are co-expressed with a specific gene of interest (in this case, 'S100A9').
```{r scSGS analysis}
SGS_res = scSGS(GetAssayData(pbmc, layer='counts'), 'S100A9', nHVG=500,
                rm.mt=T, rm.rp=T,
                calcHVG=TRUE, verbose=TRUE)

# scSGS DE results
head(SGS_res$DE)
```


# Repo contents
- /scripts - contains R files used for analysis of each dataset.
- /data - contains .xlsx containing markers genes for cell type annotation from scTypeDB.
