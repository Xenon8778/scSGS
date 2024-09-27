# single-cell Stochastic Gene Silencing (scSGS)
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
# Repo contents
- /scripts - contains R files used for analysis of each dataset.
- /data - contains .xlsx containing markers genes for cell type annotation from scTypeDB.
