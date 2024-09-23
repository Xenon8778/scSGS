futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
library(Seurat)
library(tidyverse)
library(ggpubr)
library(scSGS)

colmap_SGS = c('Silenced' = hcl.colors(4,'Peach')[2], 'Active' = 'grey')
################################
### Loading
################################
so.5K = readRDS('PBMC/PBMC5K_annot.rds')
so.10K = readRDS('PBMC/PBMC10K_annot.rds')
so.20K = readRDS('PBMC/PBMC20K_annot.rds')

# Plotting STAT1 expression
p1 = FeaturePlot(so.5K, 'STAT1', order = F)+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black'))+
  NoLegend()+ggtitle('PBMC 5K')
p2 = FeaturePlot(so.10K, 'STAT1', order = F)+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black'))+
  NoLegend()+ggtitle('PBMC 10K')
p3 = FeaturePlot(so.20K, 'STAT1', order = F)+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black'))+
  NoLegend()+ggtitle('PBMC 20K')

png('plots/PBMC/STAT1_All.png', height = 4, width = 16, units = 'in', res = 300)
ggarrange(p1,p2,p3, nrow = 1)
dev.off()

# Plotting T cells and
CTOI = so.5K$celltype
CTOI = c(rep('Other cells',length(so.5K$celltype)))
CTOI[so.5K$celltype == 'CD4 T cells'] = 'CD4 T cells'
so.5K$CTOI = factor(CTOI)

CTOI = so.10K$celltype
CTOI = c(rep('Other cells',length(so.10K$celltype)))
CTOI[so.10K$celltype == 'CD4 T cells'] = 'CD4 T cells'
so.10K$CTOI = factor(CTOI)

CTOI = so.20K$celltype
CTOI = c(rep('Other cells',length(so.20K$celltype)))
CTOI[so.20K$celltype == 'CD4 T cells'] = 'CD4 T cells'
so.20K$CTOI = factor(CTOI)

p1 = DimPlot(object = so.5K, reduction = "umap", group.by = 'CTOI', repel = T,
             label = T ,label.box = T, cols = c(hcl.colors(4,'Peach')[2], 'lightgrey'),
             label.size = 5)+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black'))+
  NoLegend()+ggtitle(NULL)
svg('plots/PBMC/UMAP_5K.svg', height = 4, width = 4)
ggarrange(p1, nrow = 1)
dev.off()
p2 = DimPlot(object = so.10K, reduction = "umap", group.by = 'CTOI', repel = T,
             label = T ,label.box = T, cols = c(hcl.colors(4,'Peach')[2], 'lightgrey'),
             label.size = 5)+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black'))+
  NoLegend()+ggtitle(NULL)
svg('plots/PBMC/UMAP_10K.svg', height = 4, width = 4)
ggarrange(p2, nrow = 1)
dev.off()
p3 = DimPlot(object = so.20K, reduction = "umap", group.by = 'CTOI', repel = T,
             label = T ,label.box = T, cols = c(hcl.colors(4,'Peach')[2], 'lightgrey'),
             label.size = 5)+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black'))+
  NoLegend()+ggtitle(NULL)
svg('plots/PBMC/UMAP_20K.svg', height = 4, width = 4)
ggarrange(p3, nrow = 1)
dev.off()

png('plots/PBMC/UMAP_All.png', height = 4, width = 12, units = 'in', res = 300)
ggarrange(p1,p2,p3, nrow = 1)
dev.off()


################
# PBMC 5K
################

## Filter MALAT1, Ribosomal and Mitochodrial genes
so_sub <- so.5K[!grepl("MALAT1", rownames(so.5K)), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'RPS'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'RPL'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'MT-'), ]
# so_sub <- so

## Subsetting
so_T = subset(so_sub, subset = celltype == c('CD4 T cells'))

## scSGS
SGS_res = scSGS(GetAssayData(so_T,layer = 'counts'),'STAT1',nHVG = 500, calcHVG = T)
so_T.SGS = so_T
so_T.SGS$GoI_SGS = SGS_res$GoI_Mask

## Differential expression
DE_wilcox_STS_Full <- SGS_res$DE %>% arrange(p_val_adj)
# feat = DE_wilcox_STS_Full %>% pull(genes)
# df = t(as.matrix(GetAssayData(so_T.SGS[feat,], layer = 'data')))
# corr_mat = cor(df,df[,'STAT1'],method = 'spearman')
# DE_wilcox_STS_Full = DE_wilcox_STS_Full %>% mutate(corr = corr_mat[feat,]) %>%
#   mutate(abs_corr = abs(corr_mat[feat,]))
head(DE_wilcox_STS_Full)

### Saving
write.csv(DE_wilcox_STS_Full,'Results/STAT1/PBMC5K_SGS.csv')
write.csv(SGS_res$HVG_df,'Results/STAT1/PBMC5K_SGS_HVG.csv')

################
# PBMC 10K
################

## Filter MALAT1, Ribosomal and Mitochodrial genes
so_sub <- so.10K[!grepl("MALAT1", rownames(so.10K)), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'RPS'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'RPL'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'MT-'), ]
# so_sub <- so

## Subsetting
so_T = subset(so_sub, subset = celltype == c('CD4 T cells'))

## scSGS
SGS_res = scSGS(GetAssayData(so_T,layer = 'counts'),'STAT1',nHVG = 500, calcHVG = T)
so_T.SGS = so_T
so_T.SGS$GoI_SGS = SGS_res$GoI_Mask

## Differential expression
DE_wilcox_STS_Full <- SGS_res$DE %>% arrange(p_val_adj)
# feat = DE_wilcox_STS_Full %>% pull(genes)
# df = t(as.matrix(GetAssayData(so_T.SGS[feat,], layer = 'data')))
# corr_mat = cor(df,df[,'STAT1'],method = 'spearman')
# DE_wilcox_STS_Full = DE_wilcox_STS_Full %>% mutate(corr = corr_mat[feat,]) %>%
#   mutate(abs_corr = abs(corr_mat[feat,]))
head(DE_wilcox_STS_Full)

### Saving
write.csv(DE_wilcox_STS_Full,'Results/STAT1/PBMC10K_SGS.csv')
write.csv(SGS_res$HVG_df,'Results/STAT1/PBMC10K_SGS_HVG.csv')

################
# PBMC 20K
################

## Filter MALAT1, Ribosomal and Mitochodrial genes
so_sub <- so.20K[!grepl("MALAT1", rownames(so.20K)), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'RPS'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'RPL'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'MT-'), ]
# so_sub <- so

## Subsetting
so_T = subset(so_sub, subset = celltype == c('CD4 T cells'))

## scSGS
SGS_res = scSGS(GetAssayData(so_T,layer = 'counts'),'STAT1',nHVG = 500, calcHVG = T)
so_T.SGS = so_T
so_T.SGS$GoI_SGS = SGS_res$GoI_Mask

## Differential expression
DE_wilcox_STS_Full <- SGS_res$DE %>% arrange(p_val_adj)
head(DE_wilcox_STS_Full)

### Saving
write.csv(DE_wilcox_STS_Full,'Results/STAT1/PBMC20K_SGS.csv')
write.csv(SGS_res$HVG_df,'Results/STAT1/PBMC20K_SGS_HVG.csv')

###############
# Comparison
###############

SGS_5K = read.csv('Results/STAT1/PBMC5K_SGS.csv')
SGS_10K = read.csv('Results/STAT1/PBMC10K_SGS.csv')
SGS_20K = read.csv('Results/STAT1/PBMC20K_SGS.csv')

sig_5K = SGS_5K %>% filter(p_val_adj< 0.01)
sig_10K = SGS_10K %>% filter(p_val_adj< 0.01)
sig_20K = SGS_20K %>% filter(p_val_adj< 0.01)

library(VennDiagram)
venn = venn.diagram(x = list('5K' = sig_5K$genes,
                             '10K' = sig_10K$genes,
                             '20K' = sig_20K$genes),
                    filename = NULL, col = 'black',
                    fill = hcl.colors(3,'peach'), cex = 1.5, cat.cex = 1.5,
                    fontfamily = "sans", ext.text = F,
                    cat.fontface = "bold",cat.fontfamily = "sans",
                    cat.dist = c(0.07,0.07,0.07))
svg('plots/PBMC/Overlaps.svg', height = 4, width = 4)
ggarrange(venn)+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "in"))
dev.off()

df = data.frame(list(dataset = factor(c('5K','10K','20K'),
                                      levels=c('5K','10K','20K')),
                     num = c(length(sig_5K$genes),
                             length(sig_10K$genes),
                             length(sig_20K$genes))))

png('plots/PBMC/bar_scSGS.png', height = 4, width = 3,units = 'in', res = 300)
ggplot(df, aes(x=dataset, y=num, fill = dataset))+
  geom_col() + theme_classic() + xlab('Datasets') +
  scale_fill_manual(values = c('steelblue2','firebrick3','goldenrod1')) +
  ylab('# scSCS responsive genes') + NoLegend()
dev.off()

# Pval and Correlation - 5K

p2 = SGS_5K %>% filter(p_val_adj < 0.05) %>% filter(genes != 'STAT1') %>%
  select(p_val_adj, abs_corr) %>% arrange(p_val_adj) %>%
  mutate(log_pval = -log(p_val_adj)) %>%
  ggplot(aes(x = abs_corr, y = log_pval))+ theme_bw() +
  geom_point(size = 1)+geom_smooth(method = lm, col = 'firebrick2')+
  stat_cor(method="pearson", label.y.npc="top", label.x.npc = "left")+
  ylab('-Log(p-value)')+xlab(expression('abs('*rho*')'))
p2
png('plots/PBMC/Correlation_5K.png', width = 4, height = 4, units = 'in', res = 300)
p2
dev.off()

# Pval and Correlation - 10K

p2 = SGS_10K %>% filter(p_val_adj < 0.05) %>% filter(genes != 'STAT1') %>%
  select(p_val_adj, abs_corr) %>% arrange(p_val_adj) %>%
  mutate(log_pval = -log(p_val_adj)) %>%
  ggplot(aes(x = abs_corr, y = log_pval))+ theme_bw() +
  geom_point(size = 1)+geom_smooth(method = lm, col = 'firebrick2')+
  stat_cor(method="pearson", label.y.npc="top", label.x.npc = "left")+
  ylab('-Log(p-value)')+xlab(expression('abs('*rho*')'))
p2
png('plots/PBMC/Correlation_10K.png', width = 4, height = 4, units = 'in', res = 300)
p2
dev.off()

# Pval and Correlation - 20K

p2 = SGS_20K %>% filter(p_val_adj < 0.05) %>% filter(genes != 'STAT1') %>%
  select(p_val_adj, abs_corr) %>% arrange(p_val_adj) %>%
  mutate(log_pval = -log(p_val_adj)) %>%
  ggplot(aes(x = abs_corr, y = log_pval))+ theme_bw() +
  geom_point(size = 1)+geom_smooth(method = lm, col = 'firebrick2')+
  stat_cor(method="pearson", label.y.npc="top", label.x.npc = "left")+
  ylab('-Log(p-value)')+xlab(expression('abs('*rho*')'))
p2
png('plots/PBMC/Correlation_20K.png', width = 4, height = 4, units = 'in', res = 300)
p2
dev.off()

###############
# Enrichr
###############
library(enrichR)
setEnrichrSite("Enrichr")
dbs <- c("GO_Biological_Process_2023")

## 5K
gl = sig_5K$genes[1:200]
SGS_Full_enrichr <- enrichr(gl, dbs)
genenames = sapply(SGS_Full_enrichr[[1]]$Genes[1:50],str_to_upper)
genenames = sapply(genenames,
                   function(x) {
                     x = gsub("\\;",", ",as.character(x))
                     gsub("^"," ",as.character(x))
                   })

Ont_5K = SGS_Full_enrichr[[1]] %>% filter(Adjusted.P.value<0.05) %>% pull(Term)

png('plots/PBMC/5K_Enrichr_GOBP.png', height = 3, width = 7, units = 'in', res = 600 )
plotEnrich(SGS_Full_enrichr[[1]], showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",
           title = 'PBMC 5K')+ scale_fill_gradient(low = '#ED90A4',high = '#5AB5E2')
dev.off()

## 10K
gl = sig_10K$genes[1:100]
SGS_Full_enrichr <- enrichr(gl, dbs)
genenames = sapply(SGS_Full_enrichr[[1]]$Genes[1:50],str_to_upper)
genenames = sapply(genenames,
                   function(x) {
                     x = gsub("\\;",", ",as.character(x))
                     gsub("^"," ",as.character(x))
                   })

Ont_10K = SGS_Full_enrichr[[1]] %>% filter(Adjusted.P.value<0.05) %>% pull(Term)

png('plots/PBMC/10K_Enrichr_GOBP.png', height = 3, width = 7, units = 'in', res = 600 )
plotEnrich(SGS_Full_enrichr[[1]], showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",
           title = 'PBMC 10K')+ scale_fill_gradient(low = '#ED90A4',high = '#5AB5E2')
dev.off()

## 20K
gl = sig_20K$genes[1:100]
SGS_Full_enrichr <- enrichr(gl, dbs)
genenames = sapply(SGS_Full_enrichr[[1]]$Genes[1:50],str_to_upper)
genenames = sapply(genenames,
                   function(x) {
                     x = gsub("\\;",", ",as.character(x))
                     gsub("^"," ",as.character(x))
                   })

Ont_20K = SGS_Full_enrichr[[1]] %>% filter(Adjusted.P.value<0.05) %>% pull(Term)

png('plots/PBMC/20K_Enrichr_GOBP.png', height = 3, width = 7, units = 'in', res = 600 )
plotEnrich(SGS_Full_enrichr[[1]], showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",
           title = 'PBMC 20K')+ scale_fill_gradient(low = '#ED90A4',high = '#5AB5E2')
dev.off()

# Intersection Human PBMCs
gl = intersect(intersect(sig_5K$genes,sig_10K$genes), sig_20K$genes)
SGS_Full_enrichr <- enrichr(gl, dbs)
genenames = sapply(SGS_Full_enrichr[[1]]$Genes[1:50],str_to_title)
genenames = sapply(genenames,
                   function(x) {
                     x = gsub("\\;",", ",as.character(x))
                     gsub("^"," ",as.character(x))
                   })

write.csv(SGS_Full_enrichr[[1]] %>% filter(P.value<0.05),'Results/STAT1/Enrichr_Intersect_SGS.csv')

## Circos plot - Chord diagram
source('SGS/Chord.R')
df = SGS_Full_enrichr[[1]]
terms = c('Response To Type II Interferon (GO:0034341)',
          'Type I Interferon-Mediated Signaling Pathway (GO:0060337)',
          'Defense Response To Virus (GO:0051607)')
make_chord(df, terms = terms, font.cex = 1, pal = 'Peach',link.trans =0.3)

# Make the circular plot
svg('plots/PBMC/SGS_STAT1_Chord.svg', height = 6, width = 6)
make_chord(df, terms = terms, font.cex = 1, pal = 'Peach', link.trans = 0.3)
dev.off()

png('plots/PBMC/Intersect_Enrichr_GOBP.png', height = 3, width = 7, units = 'in', res = 600 )
plotEnrich(SGS_Full_enrichr[[1]], showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",
           title = 'GOBP 2023')+ scale_fill_gradient(low = 'firebrick3',high = 'orange')+
  geom_text(aes(label = genenames[1:10]), colour = 'white',y = 0, size = 2, hjust = 'left',  fontface = 'bold')
dev.off()


## Pathways Compared

venn = venn.diagram(x = list('5K' = Ont_5K,
                             '10K' = Ont_10K,
                             '20K' = Ont_20K),
                             #'Mouse' = Ont_Mouse),
                    filename = NULL, col = 'black',
                    fill = hcl.colors(3,'Peach'), cex = 1.5, cat.cex = 1.5,
                    fontfamily = "sans", ext.text = F,
                    cat.fontface = "bold",cat.fontfamily = "sans",
                    cat.dist = c(0.04,0.1,0.05) )
svg('plots/PBMC/Pathway_Overlaps.svg', height = 4, width = 4)
ggarrange(venn)+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "in"))
dev.off()


##############
# Gene Violins
##############

## scSGS
int_all = intersect(intersect(sig_5K$genes,sig_10K$genes), sig_20K$genes)
clipr::write_clip((int_all))

feat = c('STAT1','STAT2','SP100','PARP14','PARP9','GBP5','GBP2',
         'GBP4','NLRC5','OAS2','IFI16','EPSTI1','EIF2AK2','RACK1',
         'APOL6','SAMD9L','XAF1','TAP1')

set.seed(1) # Selecting 6 random genes
common_genes = intersect(rownames(so.5K),intersect(rownames(so.10K),rownames(so.20K)))
featn = sample(1:length(common_genes), length(feat), replace=F)
rand_feat = common_genes[featn+1]

### 5K
SGS_res = scSGS(GetAssayData(so.5K,
                             layer = 'counts'),'STAT1',nHVG = 500)
so.5K.SGS = so.5K
so.5K.SGS$GoI_SGS = SGS_res$GoI_Mask
so.5K.SGS.T = subset(so.5K.SGS, subset = celltype == 'CD4 T cells')

p1 = VlnPlot(so.5K.SGS.T, features = feat,
             group.by = 'GoI_SGS', stack = T, flip = T, fill.by = 'ident',
             add.noise = F,cols = colmap_SGS)+ NoLegend()+
  geom_boxplot(width = 0.2, outlier.size = 0.3)+ggtitle("PBMC5K")
p1

### 10K
SGS_res = scSGS(GetAssayData(so.10K,layer = 'counts'),'STAT1',nHVG = 500)
so.10K.SGS = so.10K
so.10K.SGS$GoI_SGS = SGS_res$GoI_Mask
so.10K.SGS.T = subset(so.10K.SGS, subset = celltype == 'CD4 T cells')

p2 = VlnPlot(so.10K.SGS.T, features = feat,
             group.by = 'GoI_SGS', stack = T, flip = T, fill.by = 'ident',
             add.noise = F,cols = colmap_SGS)+NoLegend()+
  geom_boxplot(width = 0.2, outlier.size = 0.3)+ggtitle("PBMC 10K")
p2

### 20K
SGS_res = scSGS(GetAssayData(so.20K,layer = 'counts'),'STAT1',nHVG = 500)
so.20K.SGS = so.20K
so.20K.SGS$GoI_SGS = SGS_res$GoI_Mask
so.20K.SGS.T = subset(so.20K.SGS, subset = celltype == 'CD4 T cells')

p3 = VlnPlot(so.20K.SGS.T, features = feat,
             group.by = 'GoI_SGS', stack = T, flip = T, fill.by = 'ident',
             add.noise = F,cols = colmap_SGS)+NoLegend()+
  geom_boxplot(width = 0.2, outlier.size = 0.3)+ggtitle("PBMC 20K")
p3

# Dot plots
q1 = DotPlot(so.5K.SGS.T, features = feat, group.by = 'GoI_SGS',
             dot.scale = 5, cols = c('grey',hcl.colors(4,'Peach')[2]))+
  coord_flip() + NoLegend() + ylab('PBMC 5K') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
q1
q2 = DotPlot(so.10K.SGS.T, features = feat, group.by = 'GoI_SGS',
             dot.scale = 5, cols = c('grey',hcl.colors(4,'Peach')[2]))+
  coord_flip() + NoLegend() + ylab('PBMC 10K')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
q2
q3 = DotPlot(so.20K.SGS.T, features = feat, group.by = 'GoI_SGS',
             dot.scale = 5, cols = c('grey',hcl.colors(4,'Peach')[2]))+
  coord_flip() + ylab('PBMC 20K')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
q3

png('plots/PBMC/Vln_4genes.png', height = 5, width = 9, units = 'in', res = 600)
q1+q2+q3
dev.off()

# Random
p1 = VlnPlot(so.5K.SGS.T, features = rand_feat,
             group.by = 'GoI_SGS', stack = T, flip = T, fill.by = 'ident',
             add.noise = F,cols = colmap_SGS)+ NoLegend()+
  geom_boxplot(width = 0.2, outlier.size = 0.3)+ggtitle("PBMC5K")
p1

p2 = VlnPlot(so.10K.SGS.T, features = rand_feat,
             group.by = 'GoI_SGS', stack = T, flip = T, fill.by = 'ident',
             add.noise = F,cols = colmap_SGS)+NoLegend()+
  geom_boxplot(width = 0.2, outlier.size = 0.3)+ggtitle("PBMC 10K")
p2

p3 = VlnPlot(so.20K.SGS.T, features = rand_feat,
             group.by = 'GoI_SGS', stack = T, flip = T, fill.by = 'ident',
             add.noise = F,cols = colmap_SGS)+NoLegend()+
  geom_boxplot(width = 0.2, outlier.size = 0.3)+ggtitle("PBMC 20K")
p3

png('plots/PBMC/Vln_4genes_Rand.png', height = 6, width = 9, units = 'in', res = 600)
ggarrange(p1,p2,p3,nrow = 1)
dev.off()


# Random Dot plots
q1 = DotPlot(so.5K.SGS.T, features = rand_feat, group.by = 'GoI_SGS',
             dot.scale = 5, cols = c('grey',hcl.colors(4,'Peach')[2]))+
  coord_flip() + NoLegend() + ylab('PBMC 5K')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
q2 = DotPlot(so.10K.SGS.T, features = rand_feat, group.by = 'GoI_SGS',
             dot.scale = 5, cols = c('grey',hcl.colors(4,'Peach')[2]))+
  coord_flip() + NoLegend() + ylab('PBMC 10K')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
q3 = DotPlot(so.20K.SGS.T, features = rand_feat, group.by = 'GoI_SGS',
             dot.scale = 5, cols = c('grey',hcl.colors(4,'Peach')[2]))+
  coord_flip() + ylab('PBMC 20K')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

png('plots/PBMC/Vln_4genes_Rand.png', height = 5, width = 9, units = 'in', res = 600)
q1+q2+q3
dev.off()

################################
### PBMC 20k Split Test
################################
so.20K = readRDS('PBMC/PBMC20K_annot.rds')

## Filter MALAT1, Ribosomal and Mitochodrial genes
so_sub <- so.20K[!grepl("MALAT1", rownames(so.20K)), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'RPS'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'RPL'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'MT-'), ]
so_sub

## Subsetting
so_T = subset(so_sub, subset = celltype == c('CD4 T cells'))

## Random Cell Sampling
split_FC_pval <- function(seeds){
  set.seed(seeds)
  split_mask = colnames(so_T) %in% sample(colnames(so_sub), ncol(so_sub)/2)
  so_Split1 = so_T[,split_mask]
  so_Split2 = so_T[,!split_mask]

  # scSGS with STAT1
  SGS_1 = scSGS(GetAssayData(so_Split1,layer = 'counts'),'IL7R',nHVG = 500, calcHVG = F)
  SGS_2 = scSGS(GetAssayData(so_Split2,layer = 'counts'),'IL7R',nHVG = 500, calcHVG = F)

  SGS_1_DE = SGS_1$DE[-1,] %>% filter(p_val_adj < 0.05)
  SGS_2_DE = SGS_2$DE[-1,] %>% filter(p_val_adj < 0.05)

  common_genes = intersect(SGS_1_DE$genes, SGS_2_DE$genes)
  SGS_1_DE = SGS_1_DE[common_genes,]
  SGS_2_DE = SGS_2_DE[common_genes,]

  df = cbind("pval1" = -log(SGS_1_DE$p_val_adj), "pval2" = -log(SGS_2_DE$p_val_adj),
             'FC1' = SGS_1_DE$avg_log2FC, 'FC2' = SGS_2_DE$avg_log2FC)
  p1 = ggplot(df, aes(x = pval1, y = pval2)) +
    geom_point()+
    geom_smooth(formula = 'y ~ x', method = 'lm',
                color = hcl.colors(4,'Peach')[2])+
    stat_cor(cor.coef.name = 'rho',digits = 0.001,)+
    xlab('Split1') + ylab('Split2')+
    theme_classic() + ggtitle('-Log(FDR)')

  p2 = ggplot(df, aes(x = FC1, y = FC2)) +
    geom_point() +
    geom_smooth(formula = 'y ~ x', method = 'lm',
                color = hcl.colors(4,'Peach')[2]) +
    stat_cor(cor.coef.name = 'rho',digits = 0.001,) +
    xlab('Split1') + ylab('Split2')+
    theme_classic() + ggtitle('Avg Log2(FC)')
  p3 = ggarrange(p1,p2, ncol = 1, nrow = 2)
  annotate_figure(p3, top=paste('Iteration', seeds))
}

plots = lapply(1:5,split_FC_pval)

ggarrange(plotlist = plots, nrow = 1, ncol = 5)

png('plots/PBMC/split_scatters.png', height = 6, width = 15, units = 'in', res = 300)
ggarrange(plotlist = plots, nrow = 1, ncol = 5)
dev.off()
