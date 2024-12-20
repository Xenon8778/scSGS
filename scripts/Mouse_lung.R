futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
library(Seurat)
library(tidyverse)
library(ggpubr)
library(scSGS)

# Setting Colormaps
colmap_SGS = c('Silenced' = hcl.colors(4,'Peach')[2], 'Active' = 'grey')

############
# scSGS ####
############
# Loading Data
so.raw <- readRDS('Twocondition/Lung_cancer_UMAP.rds')

# Find marker genes for cell types
Idents(so.raw) <- 'celltype'
markers <- FindAllMarkers(so.raw, only.pos = TRUE, logfc.threshold = 1)
top20 <- markers %>%
  group_by(cluster) %>%
  slice_head(n = 20) %>%
  ungroup()
write.csv(top20, 'results/Twocondition/top20markers.csv')
FeaturePlot(so.raw, 'Aplnr')

so <- subset(so.raw, subset = celltype == 'gCap') # subset gCap only
so <- RunUMAP(so, dims = 1:40)
colnames(so@meta.data)

CTOI = c(rep('gCap',length(so$celltype)))
names(CTOI) = names(so$celltype)
so.raw$CTOI = CTOI
so.raw$CTOI[is.na(so.raw$CTOI)] = 'Other cells'
so.raw$CTOI = factor(so.raw$CTOI)

p2 = DimPlot(so.raw, reduction = "umap", group.by = 'CTOI', label = T, label.box = T,
             cols = c(hcl.colors(4,'peach')[2],'lightgrey'),label.size = 5)+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black'))+
  NoLegend()+ggtitle(NULL)
p2
svg('plots/Twocondition/UMAP_Lung_CTOI.svg', height = 4, width = 4)
p2
dev.off()

p1 = DimPlot(so, reduction = "umap", group.by = 'sample', label = T,
        label.box = T, cols = c(hcl.colors(4,'Peach')[2], 'grey')) + NoLegend()+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black'))+
  NoLegend()+ggtitle(NULL)
p1
svg('plots/Twocondition/UMAP_Lung_sample.svg', height = 4, width = 4)
p1
dev.off()
p1 = FeaturePlot(so.raw, c('Plvap','Car4','Amigo2','Col3a1'), alpha = 0.3,
                 max.cutoff = 3, min.cutoff = 1, ncol = 2)
p1
FeaturePlot(so, features = c('Stk11'))

# Splintting samples
Control = subset(so, subset = sample == 'Control')
Cancer = subset(so, subset = sample == 'Cancer')

HVG_Control = HVG_splinefit(GetAssayData(Control, layer = 'counts'),
                            dropout.filter = T, nHVGs = 2000)
HVG_Cancer = HVG_splinefit(GetAssayData(Cancer, layer = 'counts'),
                            dropout.filter = T, nHVGs = 2000)
write.csv(HVG_Control,'Results/Twocondition/Control_SGS_HVG.csv')
write.csv(HVG_Cancer,'Results/Twocondition/Cancer_SGS_HVG.csv')

p1 <- FeaturePlot(Control, 'Stk11', order = F,
                 pt.size = 0.75, cols = c('lightgrey',hcl.colors(4,'Peach')[1]))+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black',linewidth = 1))+
  NoLegend() + ggtitle(NULL)
p1
p2 <- FeaturePlot(Cancer, 'Stk11', order = F,
                 pt.size = 0.75, cols = c('lightgrey',hcl.colors(4,'Peach')[1]))+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black',linewidth = 1))+
  NoLegend() + ggtitle(NULL)
p2
png('plots/Twocondition/Stk11_gCap.png', height = 6, width = 3, units = 'in', res = 300)
ggarrange(p1,p2, ncol = 1)
dev.off()

Common_HVG = intersect(rownames(HVG_Control[HVG_Control$HVG == T,]),
                       rownames(HVG_Cancer[HVG_Cancer$HVG == T,]))
head(Common_HVG)
Common_HVG[startsWith(Common_HVG,'Stk')]

## scSGS Control ####
SGS_Control = scSGS(GetAssayData(Control,layer = 'counts'),'Stk11',nHVG = 500,
                    rm.mt = T, rm.rp = T, calcHVG = F)
SGS_Control.SGS = Control
SGS_Control.SGS$GoI_SGS = SGS_Control$GoI_Mask
DE_SGS_Control <- SGS_Control$DE %>% arrange(p_val_adj)

DimPlot(SGS_Control.SGS, reduction = "umap", group.by = 'GoI_SGS')+
  ggtitle('Control')

## scSGS Cancer ####
SGS_Cancer = scSGS(GetAssayData(Cancer,layer = 'counts'),'Stk11',nHVG = 500,
                   rm.mt = T, rm.rp = T, calcHVG = F)
SGS_Cancer.SGS = Cancer
SGS_Cancer.SGS$GoI_SGS = SGS_Cancer$GoI_Mask
DE_SGS_Cancer <- SGS_Cancer$DE %>% arrange(p_val_adj)
DimPlot(SGS_Cancer.SGS, reduction = "umap", group.by = 'GoI_SGS')+
  ggtitle('Cancer')

write.csv(DE_SGS_Control,'Results/Twocondition/SGS_Control.csv')
write.csv(DE_SGS_Cancer,'Results/Twocondition/SGS_Cancer.csv')

##################
# Comparisons ####
##################
library(VennDiagram)
library(ggplot2)
library(ggpubr)

DE_SGS_Control = read.csv('Results/Twocondition/SGS_Control.csv')
DE_SGS_Cancer = read.csv('Results/Twocondition/SGS_Cancer.csv')

nSig_Control = DE_SGS_Control %>% subset(p_val_adj<0.01) %>%
  subset(abs(avg_log2FC) > 0.25)
t200_Control = nSig_Control %>% head(200)
clipr::write_clip(t200_Control$genes[1:30])

nSig_Cancer = DE_SGS_Cancer %>% subset(p_val_adj<0.01)%>%
  subset(abs(avg_log2FC) > 0.25)
t200_Cancer = nSig_Cancer %>% head(200)
clipr::write_clip(t200_Cancer$genes)

venn1 = venn.diagram(x = list('Control' = nSig_Control$genes,
                              'Cancer' = nSig_Cancer$genes),
                     filename = NULL, col = 'black', alpha = 0.7,
                     fill = c('white',hcl.colors(4,'Peach')[2]), cex = 1.5, cat.cex = 1.5,
                     fontfamily = "sans",ext.text = F,ext.dist = -0.1, ext.length = 0.9,
                     cat.fontface = "bold",cat.fontfamily = "sans",
                     cat.dist = c(0.09,0.08), margin = 0.1)
ggarrange(venn1)
svg('plots/Twocondition/Venn_Lung_Stk11.svg', height = 4, width = 4)
ggarrange(venn1)
dev.off()

#############
# Enrichr
#############
library(enrichR)
setEnrichrSite("Enrichr")
dbs <- c("GO_Biological_Process_2023")

## Control
gl = t200_Control$genes
Control_enrichr <- enrichr(gl, dbs)
Ont_Control = Control_enrichr[[1]] %>% filter(Adjusted.P.value<0.05) %>%
  mutate(totalgenes = as.numeric(lapply(strsplit(Overlap,'/'),'[[',2))) %>%
  filter(totalgenes < 500) %>%
  mutate(ngenes = lengths(strsplit(Genes, ';')))
write.csv(Ont_Control, 'Results/Twocondition/SGS_Control_Enrichr.csv')

p1 = Ont_Control %>% ggplot(aes(x = log(Odds.Ratio), y = reorder(Term, -Adjusted.P.value),
                          size = ngenes, fill = Adjusted.P.value))+
  geom_point(stroke = 0.5, color = 'black', shape=21)+
  theme_classic()+
  labs(y = 'GOBP', x = 'log(OR)', fill = 'adj.P', size = 'Count')+
  scale_fill_gradientn(colours = hcl.colors(9, 'peach'))+
  lims(size = c(1,NA))+
  ggtitle('Enriched Pathways')
p1

png('plots/Twocondition/Bubble_Lung_Control_Stk11.png', height = 4, width = 8.5, res = 300, units = 'in')
p1
dev.off()

## Cancer
gl = t200_Cancer$genes
Cancer_enrichr <- enrichr(gl, dbs)
Ont_Cancer = Cancer_enrichr[[1]] %>% filter(Adjusted.P.value<0.05) %>%
  mutate(totalgenes = as.numeric(lapply(strsplit(Overlap,'/'),'[[',2))) %>%
  filter(totalgenes < 500) %>%
  mutate(ngenes = lengths(strsplit(Genes, ';')))
write.csv(Ont_Cancer, 'Results/Twocondition/SGS_Cancer_Enrichr.csv')

p2 = Ont_Cancer %>% head(15) %>%
  ggplot(aes(x = log(Odds.Ratio), y = reorder(Term, -Adjusted.P.value),
             size = ngenes, fill = Adjusted.P.value))+
  geom_point(stroke = 0.5, color = 'black', shape=21)+
  theme_classic()+
  labs(y = 'GOBP', x = 'log(OR)', fill = 'adj.P', size = 'Count')+
  scale_fill_gradientn(colours = hcl.colors(9, 'peach'))+
  lims(size = c(1,14))+
  ggtitle('Enriched Pathways')
p2

png('plots/Twocondition/Bubble_Lung_Cancer_Stk11.png', height = 4, width = 8, res = 300, units = 'in')
p2
dev.off()

Uniq_path <- setdiff(Ont_Cancer$Term, Ont_Control$Term)
Ont_Uniq <- Ont_Cancer %>% filter(Term %in% Uniq_path)

p3 = Ont_Uniq %>% head(15) %>%
  ggplot(aes(x = log(Odds.Ratio), y = reorder(Term, -Adjusted.P.value),
             size = ngenes, fill = Adjusted.P.value))+
  geom_point(stroke = 0.5, color = 'black', shape=21)+
  theme_classic()+
  labs(y = 'GOBP', x = 'log(OR)', fill = 'adj.P', size = 'Count')+
  scale_fill_gradientn(colours = hcl.colors(9, 'peach'))+
  lims(size = c(0,8))+
  ggtitle('Enriched Pathways')
p3

png('plots/Twocondition/Bubble_Lung_Uniq_Stk11.png', height = 4, width = 7.5, res = 300, units = 'in')
p3
dev.off()

# Pathway overlap
venn1 = venn.diagram(x = list('Control' = Ont_Control$Term,
                              'Cancer' = Ont_Cancer$Term),
                     filename = NULL, col = 'black', alpha = 0.7,
                     fill = c('white',hcl.colors(4,'Peach')[2]), cex = 1.5, cat.cex = 1.5,
                     fontfamily = "sans",ext.text = F,ext.dist = -0.1, ext.length = 0.9,
                     cat.fontface = "bold",cat.fontfamily = "sans",
                     cat.dist = c(0.08,0.07), margin = 0.1)
ggarrange(venn1)

svg('plots/Twocondition/Venn_Lung_Pathway_Stk11.svg', height = 4, width = 4)
ggarrange(venn1)
dev.off()

## SGS Uniq Genes
gl = setdiff(t200_Cancer$genes, t200_Control$genes) %>% head(200)
Uniq_enrichr <- enrichr(gl, dbs)
p3 = plotEnrich(Uniq_enrichr[[1]], numChar = 80)
p3
Ont_Uniq = Uniq_enrichr[[1]] %>% filter(Adjusted.P.value<0.05)

p1 + p2

## Chords
source('SGS/Chord.R')

df = Ont_Uniq
terms = c('Anoikis (GO:0043276)',
          'Regulation Of NIK/NF-kappaB Signaling (GO:1901222)',
          'Regulation Of Signal Transduction By P53 Class Mediator (GO:1901796)',
          'Vasculogenesis (GO:0001570)')
make_chord(df, terms = terms, font.cex = 1, pal = 'Peach',link.trans =0.3,
           small.gap = 0.2)





