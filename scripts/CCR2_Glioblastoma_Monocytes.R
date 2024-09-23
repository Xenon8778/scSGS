futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
library(Seurat)
library(tidyverse)
library(scSGS)

###################
## RNA
###################
# Loading Data
so = readRDS('Ccr2/RNA/RNA_glioblastoma_UMAP.rds')
so$Batches = so$sample %>% substr(start = 1 ,stop = 2) # Merge replicates

## Filter MALAT1, Ribosomal and Mitochodrial genes
so_sub <- so[!grepl("Malat1", rownames(so)), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'Rps'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'Rpl'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'mt-'), ]

## Subsetting
#so_TAM1 = subset(so_sub, subset = cluster == c('TAM 1'))
so_Monocytes = subset(so_sub, subset = cluster == c('Monocytes')) # Extract monocytes
so_Whole = so_sub

###################
## Whole
###################
WT_whole = subset(so_Whole, subset = Batches == "WT")
# scSGS analysis
SGS_res = scSGS(GetAssayData(WT_whole,layer = 'counts'),'Ccr2',nHVG = 500)
so_WT.SGS = WT_whole
so_WT.SGS$GoI_SGS = SGS_res$GoI_Mask
# Cellcycle scoring
so_WT.SGS = CellCycleScoring(so_WT.SGS, s.features = cc.genes$s.genes,
                             g2m.features = cc.genes$g2m.genes)
table(so_WT.SGS$Phase)

## selecting only G1 and S phase
so_G1 = subset(so_WT.SGS, subset = Phase == 'G1') # Extract G1 phase cells
SGS_res_G1 = scSGS(GetAssayData(so_G1,layer = 'counts'),'Ccr2') # scSGS analysis
so_S = subset(so_WT.SGS, subset = Phase == 'S') # Extract S phase cells
SGS_res_S = scSGS(GetAssayData(so_S,layer = 'counts'),'Ccr2') # scSGS analysis

## Differential expression
### Full
DE_wilcox_STS_Full <- get.DE(SGS_res) %>% arrange(p_val_adj)
feat = DE_wilcox_STS_Full %>% pull(genes)
df = t(as.matrix(GetAssayData(so_WT.SGS[feat,], layer = 'data')))
corr_mat = cor(df,df[,'Ccr2'],method = 'spearman')
DE_wilcox_STS_Full = DE_wilcox_STS_Full %>% mutate(corr = corr_mat[feat,]) %>%
  mutate(abs_corr = abs(corr_mat[feat,]))
head(DE_wilcox_STS_Full)

### G1 Phase
DE_wilcox_G1 <- get.DE(SGS_res_G1) %>% arrange(p_val_adj)
feat = DE_wilcox_G1 %>% pull(genes)
df = t(as.matrix(GetAssayData(so_G1[feat,], layer = 'data')))
corr_mat = cor(df,df[,'Ccr2'],method = 'spearman')
DE_wilcox_G1 = DE_wilcox_G1 %>% mutate(corr = corr_mat[feat,]) %>%
  mutate(abs_corr = abs(corr_mat[feat,]))
head(DE_wilcox_G1)

### S Phase
DE_wilcox_S <- get.DE(SGS_res_S) %>% arrange(p_val_adj)
feat = DE_wilcox_S %>% pull(genes)
df = t(as.matrix(GetAssayData(so_S[feat,], layer = 'data')))
corr_mat = cor(df,df[,'Ccr2'],method = 'spearman')
DE_wilcox_S = DE_wilcox_S %>% mutate(corr = corr_mat[feat,]) %>%
  mutate(abs_corr = abs(corr_mat[feat,]))
head(DE_wilcox_S)

### Saving
write.csv(DE_wilcox_G1,'Results/Ccr2/Whole/Whole_SGS_DE_G1.csv')
write.csv(DE_wilcox_S,'Results/Ccr2/Whole/Whole_SGS_DE_S.csv')
write.csv(DE_wilcox_STS_Full,'Results/Ccr2/Whole/Whole_SGS_DE_Full.csv')

## In vivo KO
### Pre-processing
so_Whole <- NormalizeData(so_Whole)
res = FindMarkers(so_Whole, ident.1 = 'WT', ident.2 = 'KO', group.by = 'Batches',
                  logfc.threshold = 0.1, min.pct = 0.1, test.use = "wilcox")
### DE
DE_wilcox_OG <- res %>% mutate(genes = row.names(res)) %>% arrange(p_val_adj)
head(DE_wilcox_OG)
write.csv(DE_wilcox_OG,'Results/Ccr2/Whole/Whole_SGS_DE_OG.csv')

###################
## Monocytes only
###################

WT_Mon = subset(so_Monocytes, subset = Batches == "WT")
# scSGS analysis
SGS_res = scSGS(GetAssayData(WT_Mon,layer = 'counts'),'Ccr2',nHVG = 500, calcHVG = T)
WT_Mon.SGS = WT_Mon
WT_Mon.SGS$GoI_SGS = SGS_res$GoI_Mask
WT_Mon.SGS = CellCycleScoring(WT_Mon.SGS, s.features = cc.genes$s.genes,
                             g2m.features = cc.genes$g2m.genes)
table(WT_Mon.SGS$Phase)

## selecting only G1 phase
so_G1 = subset(WT_Mon.SGS, subset = Phase == 'G1') # Extract G1 phase cells
SGS_res_G1 = scSGS(GetAssayData(so_G1,layer = 'counts'),'Ccr2')# scSGS analysis
so_S = subset(WT_Mon.SGS, subset = Phase == 'S') # Extract S phase cells
SGS_res_S = scSGS(GetAssayData(so_S,layer = 'counts'),'Ccr2')# scSGS analysis


## Differential expression
### Full
DE_wilcox_STS_Full <- SGS_res$DE %>% arrange(p_val_adj)
feat = DE_wilcox_STS_Full %>% pull(genes)
df = t(as.matrix(GetAssayData(WT_Mon.SGS[feat,], layer = 'data')))
head(DE_wilcox_STS_Full)


### G1 Phase
DE_wilcox_G1 <- SGS_res_G1$DE %>% arrange(p_val_adj)
feat = DE_wilcox_G1 %>% pull(genes)
df = t(as.matrix(GetAssayData(so_G1[feat,], layer = 'data')))
head(DE_wilcox_G1)


### S Phase
DE_wilcox_S <- SGS_res_S$DE %>% arrange(p_val_adj)
feat = DE_wilcox_S %>% pull(genes)
df = t(as.matrix(GetAssayData(so_S[feat,], layer = 'data')))
head(DE_wilcox_S)

### Saving
write.csv(DE_wilcox_G1,'Results/Ccr2/Monocytes/Monocytes_SGS_DE_G1.csv')
write.csv(DE_wilcox_S,'Results/Ccr2/Monocytes/Monocytes_SGS_DE_S.csv')
write.csv(DE_wilcox_STS_Full,'Results/Ccr2/Monocytes/Monocytes_SGS_DE_Full.csv')

write.csv(SGS_res$HVG_df,'Results/Ccr2/Monocytes/Monocytes_SGS_HVG_Full.csv')

## Invivo KO
### Pre-processing
so_Monocytes <- NormalizeData(so_Monocytes)
### DE
res = FindMarkers(so_Monocytes, ident.1 = 'KO', ident.2 = 'WT', group.by = 'Batches',
                  logfc.threshold = 0, min.pct = 0.1, test.use = "wilcox")
DE_wilcox_OG <- res %>% mutate(genes = row.names(res)) %>% arrange(p_val_adj)
head(DE_wilcox_OG)
write.csv(DE_wilcox_OG,'Results/Ccr2/Monocytes/Monocytes_SGS_DE_OG.csv')

