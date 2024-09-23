# single-cell Stochastic Gene Silencing (scSGS)
![Figure 1 - Compressed](https://github.com/Xenon8778/scSGS/assets/61325229/9c629e1e-4b34-456d-a80e-7476622ff6b4)

This is the main repository for scSGS - A method to study gene function using the natural stochasticity in gene expression.

The scripts for preprocessing and analyzing publicly available datasets are in the '/scripts' folder.

# Install scSGS
```R
# install.packages("devtools")
devtools::install_github("immunogenomics/presto")
devtools::install_github("Xenon8778/scSGS")
```
# Repo contents
"/scripts" - contains R files used for analysis of each dataset.
"/data" - contains .xlsx containing markers genes for cell type annotation from scTypeDB.
