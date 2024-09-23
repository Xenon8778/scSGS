futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
library(VennDiagram)
library(ggpubr)
library(Seurat)
library(tidyverse)
library(scSGS)
library(ggrepel)

# PBMC 20K ####
so = readRDS('PBMC/PBMC20K_annot.rds')

## Filter MALAT1, Ribosomal and Mitochodrial genes
so_sub <- so[!grepl("MALAT1", rownames(so)), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'RPS'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'RPL'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'MT-'), ]

## Extract CD4 T cells and top 5000 HVG
so_T = subset(so_sub, subset = celltype == c('CD4 T cells'))
HVGs = HVG_splinefit(GetAssayData(so_T, layer = 'count'), nHVGs = 5000)
write.csv(HVGs, 'Results/Correlation/HVG_IL7R_PBMC20K.csv')

hvgenes = rownames(HVGs %>% filter(HVG == T))
so_T_5000 <- so_T[hvgenes,]

## scSGS analysis ####
scSGS_res <- scSGS(GetAssayData(so_T_5000, layer = 'count'), GoI = 'IL7R')
SGS_DF <- scSGS_res$DE
write.csv(SGS_DF, 'Results/Correlation/SGS_IL7R_PBMC20K.csv')

####################
# Spearman's Correlation
####################
so_T_mat = t(as.matrix(GetAssayData(so_T_5000, layer = 'count')))
cor_spearman = cor(so_T_mat[,'IL7R'],so_T_mat, method = 'spearman')
cor_spearman_df = as.data.frame(t(cor_spearman))

# Compute P-values from Spearman
cor_spearman_pval = apply(so_T_mat, 2, function(x) {
  cor.test(so_T_mat[,'IL7R'], x, method = "spearman")$p.value
})
cor_spearman_df_pval = as.data.frame(cor_spearman_pval)

SGS_DF['spearman'] = cor_spearman_df[row.names(SGS_DF),]
SGS_DF['spearman_pval'] = cor_spearman_df_pval[row.names(SGS_DF),]
SGS_DF['spearman_FDR'] = p.adjust(SGS_DF$spearman_pval)

p1 = SGS_DF[-1,] %>% ggplot() +
  geom_point(aes(x = -avg_log2FC, y = spearman, color = p_val_adj < 0.01),
             shape = 16, alpha = 0.7, size = 0.7)+
  geom_point(SGS_DF[-1,] %>% filter(p_val_adj < 0.01),
             mapping = aes(x = -avg_log2FC, y = spearman),
             shape = 16, size = 1.5, color = 'firebrick3')+
  geom_point(SGS_DF[-1,] %>% filter(p_val_adj < 0.01 & spearman_FDR > 0.01),
             mapping = aes(x = -avg_log2FC, y = spearman),
             shape = 16, size = 2, color = 'orange')+
  annotate('text',label = paste0('Significant\nGenes =\n',
                                 nrow(SGS_DF %>% filter(p_val_adj < 0.01))),
           x = -3.5, y = 0.3,hjust = 0,size = 4)+
  geom_label_repel(SGS_DF[-1,] %>% filter(p_val_adj < 0.01 & spearman_FDR > 0.01) %>%
               arrange(p_val_adj) %>% head(10),
             mapping = aes(x = -avg_log2FC, y = spearman, label = genes),
             arrow = arrow(length = unit(0.2,'cm')),force = 100,size = 3)+
  scale_color_manual(values = c('black','orange'))+
  theme_classic()+ylab("Spearman's Coeff")+
  xlab('-AvgLog2FC')+NoLegend()
p1
svg('plots/Correlation/Spearman_scSGS_genes.svg', height = 5, width = 5)
p1
dev.off()

p2 = SGS_DF[-1,] %>% ggplot() +
  geom_point(aes(x = -avg_log2FC, y = spearman, color = spearman_FDR < 0.01),
             shape = 16, alpha = 0.7, size = 0.7)+
  geom_point(SGS_DF[-1,] %>% filter(spearman_FDR < 0.01),
             mapping = aes(x = -avg_log2FC, y = spearman),
             shape = 16, size = 1, color = 'firebrick3')+
  annotate('text',label = paste0('Significant\nGenes =\n',
                                 nrow(SGS_DF[-1,] %>% filter(spearman_FDR < 0.01))),
           x = -3.5, y = 0.3,hjust = 0)+
  theme_classic()+ylab("Spearman's Coeff")+
  scale_color_manual(values = c('black','firebrick3'))+
  xlab('-AvgLog2FC')+ labs(color = 'FDR<0.01')+
  NoLegend()
p1+p2

svg('plots/Correlation/Spearman_signif.svg', height = 3, width = 3)
p2
dev.off()

####################
# Pearson's Correlation
####################
so_T_mat = t(as.matrix(GetAssayData(so_T_5000, layer = 'count')))
cor_pearson = cor(so_T_mat[,'IL7R'],so_T_mat, method = 'pearson')
cor_pearson_df = as.data.frame(t(cor_pearson))

# Compute P-values from pearson
cor_pearson_pval = apply(so_T_mat, 2, function(x) {
  cor.test(so_T_mat[,'IL7R'], x, method = "pearson")$p.value
})
cor_pearson_pval = as.data.frame(cor_pearson_pval)

SGS_DF['pearson'] = cor_pearson_df[row.names(SGS_DF),]
SGS_DF['pearson_pval'] = cor_pearson_pval[row.names(SGS_DF),]
SGS_DF['pearson_FDR'] = p.adjust(SGS_DF$pearson_pval)

p1 = SGS_DF[-1,] %>% ggplot() +
  geom_point(aes(x = -avg_log2FC, y = pearson, color = p_val_adj < 0.01),
             shape = 16, alpha = 0.7, size = 0.7)+
  geom_point(SGS_DF[-1,] %>% filter(p_val_adj < 0.01),
             mapping = aes(x = -avg_log2FC, y = pearson),
             shape = 16, size = 1.5, color = 'firebrick3')+
  geom_point(SGS_DF[-1,] %>% filter(p_val_adj < 0.01 & pearson_FDR > 0.01),
             mapping = aes(x = -avg_log2FC, y = pearson),
             shape = 16, size = 2, color = 'orange')+
  annotate('text',label = paste0('Significant\nGenes =\n',
                                 nrow(SGS_DF %>% filter(p_val_adj < 0.01))),
           x = -3.5, y = 0.3,hjust = 0,size = 4)+
  geom_label_repel(SGS_DF[-1,] %>% filter(p_val_adj < 0.01 & pearson_FDR > 0.01) %>%
                     arrange(p_val_adj) %>% head(10),
                   mapping = aes(x = -avg_log2FC, y = pearson, label = genes),
                   arrow = arrow(length = unit(0.2,'cm')),force = 100,size = 3)+
  scale_color_manual(values = c('black','orange'))+
  theme_classic()+ylab("Pearson's Coeff")+
  xlab('-AvgLog2FC')+NoLegend()
p1
svg('plots/Correlation/Pearson_scSGS_genes.svg', height = 5, width = 5)
p1
dev.off()
p2 = SGS_DF[-1,] %>% ggplot() +
  geom_point(aes(x = -avg_log2FC, y = pearson, color = pearson_FDR < 0.01),
             shape = 16, alpha = 0.7, size = 0.7)+
  geom_point(SGS_DF[-1,] %>% filter(pearson_FDR < 0.01),
             mapping = aes(x = -avg_log2FC, y = pearson),
             shape = 16, size = 1, color = 'firebrick3')+
  annotate('text',label = paste0('Significant\nGenes =\n',
                                 nrow(SGS_DF[-1,] %>% filter(pearson_FDR < 0.01))),
           x = -3.5, y = 0.27,hjust = 0)+
  theme_classic()+ylab("Pearson's Coeff")+
  scale_color_manual(values = c('black','firebrick3'))+
  xlab('-AvgLog2FC')+ labs(color = 'Pvalue<0.01')+
  NoLegend()
p1+p2

svg('plots/Correlation/pearson_signif.svg', height = 3, width = 3)
p2
dev.off()

############
# Gene Set Enrichment
############
# Enrichr
top100_SGS <- SGS_DF %>% arrange(p_val_adj) %>% head(100) %>% pull(genes)
top100_spearman <- SGS_DF %>% arrange(spearman_FDR) %>% head(100) %>% pull(genes)
top100_pearson <- SGS_DF %>% arrange(pearson_FDR) %>% head(100) %>% pull(genes)

library(enrichR)
setEnrichrSite('EnrichR')
dbs = c('GO_Biological_Process_2023')
SGS_uniq = setdiff(top100_SGS, union(top100_spearman, top100_pearson))
enrichr_SGS = enrichr(SGS_uniq, dbs)[[1]]

# SGS
enrichr_SGS = enrichr(top100_SGS, dbs)[[1]] %>% filter(Adjusted.P.value < 0.05)
clipr::write_clip(enrichr_SGS)
write.csv(enrichr_SGS, 'Results/Correlation/SGS_GOBP_enrichr.csv')

# Spearman
enrichr_spearman = enrichr(top100_spearman, dbs)[[1]] %>% filter(Adjusted.P.value < 0.05)

# Pearson
enrichr_pearson = enrichr(top100_pearson, dbs)[[1]] %>% filter(Adjusted.P.value < 0.05)

venn1=venn.diagram(x=list(Spearman=enrichr_spearman$Term,
                          Pearson=enrichr_pearson$Term,
                          SGS=enrichr_SGS$Term),
                   filename=NULL, col='black', alpha=0.7,
                   fill=c('white','white',hcl.colors(4,'Peach')[2]),
                   fontfamily="sans", cat.dist=c(0.18,0.20,0.1), cex=1.5, cat.cex=1.5,
                   cat.fontface="bold",cat.fontfamily="sans",
                   margin = 0.15)
ggarrange(venn1)
svg('plots/Correlation/Venn_pearson_spearman_enrichr.svg', height=4, width=4)
ggarrange(venn1)
dev.off()

plotEnrich(enrichr_SGS,y = 'Count')

SGS_term_Unique = setdiff(enrichr_SGS$Term,
                          union(enrichr_spearman$Term, enrichr_pearson$Term))


svg('plots/Correlation/SGS_Unique_Terms.svg', height=2.5, width=5)
enrichr_SGS %>% filter(Term %in% SGS_term_Unique) %>% ggplot()+
  geom_col(aes(x = -log(Adjusted.P.value), y = reorder(Term,-Adjusted.P.value),
                fill = -log(Adjusted.P.value)), width = 0.7)+
  theme_classic()+ xlab('-log(FDR)')+ ylab('GO Terms')+
  scale_fill_gradientn(colours = hcl.colors(4,'Peach', rev = T))+
  theme(axis.text.x = element_blank())+
  NoLegend()
dev.off()

clipr::write_clip(enrichr_SGS %>% filter(Term %in% SGS_term_Unique))

