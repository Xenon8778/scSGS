# single-cell Stochastic Gene Silencing (scSGS)
Read the full article [**Here**](https://doi.org/10.1038/s42003-025-07530-0)

![Figure 1 - Compressed](https://github.com/Xenon8778/scSGS/assets/61325229/9c629e1e-4b34-456d-a80e-7476622ff6b4)

This is the main repository for scSGS - A method to study gene function using the natural stochasticity in gene expression.

The scripts for preprocessing and analyzing publicly available datasets are in the '/scripts' folder.

# Install scSGS
scSGS requires [presto](https://github.com/immunogenomics/presto) for fast computation of DE genes using the wilcoxon rank-sum test.
```R
if(!require(remotes)) install.packages('remotes')
remotes::install_github("immunogenomics/presto")
remotes::install_github("Xenon8778/scSGS")
```
# Input
scSGS needs a scRNA-seq Matrix or dgCMatrix object with genes as rows and cells as columns.

# Example
## Identification of HVGs
Identify HVGs using the `HVG_splinefit` function. This method leverages a spline-based approach to account for non-linear relationships between mean expression, coefficient of variation, and dropout rate.
```R
HVG_res = HVG_splinefit(GetAssayData(adata, layer='counts'), diptest=TRUE,
                        dropout.filter=TRUE, nHVGs=500)
```                        
## scSGS analysis
scSGS analysis should be performed for highly variable genes only. The input for the algorithm is a  scRNA-seq geneexpression matrix with genes as rows and cells as columns. The output of the algorithm is data.frame find the SGS-responsive genes and their P-values. HVG selection is done using HVG_splinefit which uses each gene's mean expression, CV, and dropout rate to compute variability.
```R
SGS_res <- scSGS(GetAssayData(adata, layer='counts'), GoI = 'LYZ', calcHVG = T)
```
## Citation
Gupta, S., Cai, J.J. Gene function revealed at the moment of sitochastic gene silencing. Commun Biol 8, 88 (2025). https://doi.org/10.1038/s42003-025-07530-0

# Repo contents
- /scripts - contains R files used for analysis of each dataset.
- /data - contains .xlsx containing markers genes for cell type annotation from scTypeDB.
