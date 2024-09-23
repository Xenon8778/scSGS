futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
library(Seurat)
library(tidyverse)
library(ggpubr)
library(scSGS)

# Setting Colormaps
colmap_Batches = c('KO' = hcl.colors(4,'Peach')[2], 'WT' = 'grey')
colmap_SGS = c('Silenced' = hcl.colors(4,'Peach')[2], 'Active' = 'grey')

# Loading Data
so = readRDS('Kdm6b/Kdm6b_Neuron.rds')
DimPlot(so, reduction = "umap", group.by = 'celltype', label = T,
        label.box = T) + NoLegend()

# identifying Slc18a3+ Clusters
so_sub = subset(so, Slc18a3 > 0)
so_sub = RunUMAP(so_sub, dims = 1:40)
p1 = FeaturePlot(so_sub, 'Kdm6b', order = T)
p2 = DimPlot(object = so_sub, reduction = "umap", group.by = 'Batch')
p1+p2

table(so_sub$Batch)
write.csv(table(so$Batch,so$celltype),'Kdm6b/CT_numbers.csv')

CTOI = c(rep('MN Lineage',length(so_sub$celltype)))
names(CTOI) = names(so_sub$celltype)
so$CTOI = CTOI
so$CTOI[is.na(so$CTOI)] = 'Other cells'
so$CTOI = factor(so$CTOI)

p2 = DimPlot(so, reduction = "umap", group.by = 'CTOI', label = T, label.box = T,
             cols = c(hcl.colors(4,'peach')[2],'lightgrey'),label.size = 5)+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black'))+
  NoLegend()+ggtitle(NULL)
p2
svg('plots/Kdm6b/UMAP_Neuron.svg', height = 4, width = 4)
p2
dev.off()

p1 = FeaturePlot(subset(so_sub, subset = Batch == 'WT'), 'Kdm6b', order = T,
                 pt.size = 1, cols = c('lightgrey',hcl.colors(4,'Peach')[1]))+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black',linewidth = 1))+
  NoLegend() + ggtitle(NULL)
p1
svg('plots/Kdm6b/Kdm6b_Neuron_WT.svg', height = 4, width = 4)
p1
dev.off()


p2 = DimPlot(object = so, reduction = "umap", group.by = 'celltype', repel = T,
             label = T ,label.box = T, cols = hcl.colors(13,'Set2'),label.size = 3)+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black'))+
  NoLegend()+ggtitle(NULL)
p2
p3 = DimPlot(object = so_sub, reduction = "umap", group.by = 'Batch', repel = T,
             label = T ,label.box = T, cols = colmap_Batches,label.size = 5)+
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(color = 'black'))+
  NoLegend()+ggtitle(NULL)
p3

## Filter MALAT1, Ribosomal and Mitochodrial genes
so_sub <- so_sub[!grepl("Malat1", rownames(so_sub)), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'Rps'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'Rpl'), ]
so_sub <- so_sub[!startsWith(rownames(so_sub),'mt-'), ]
# so_sub <- so

## Subsetting
so_Whole = so_sub

## scSG Whole ####

WT_whole = subset(so_Whole, subset = Batch == "WT")
SGS_res = scSGS(GetAssayData(WT_whole,layer = 'counts'),'Kdm6b',nHVG = 500, calcHVG = T)
so_WT.SGS = WT_whole
so_WT.SGS$GoI_SGS = SGS_res$GoI_Mask
so_WT.SGS = CellCycleScoring(so_WT.SGS, s.features = cc.genes$s.genes,
                             g2m.features = cc.genes$g2m.genes)

table(so_WT.SGS$Phase)

## selecting only G1 phase
so_G1 = subset(so_WT.SGS, subset = Phase == 'G1')
SGS_res_G1 = scSGS(GetAssayData(so_G1,layer = 'counts'),'Kdm6b',nHVG = 500)

so_S = subset(so_WT.SGS, subset = Phase == 'S')
SGS_res_S = scSGS(GetAssayData(so_S,layer = 'counts'),'Kdm6b',nHVG = 500)

## Differential expression
### Full
DE_wilcox_STS_Full <- SGS_res$DE %>% arrange(p_val_adj)
feat = DE_wilcox_STS_Full %>% pull(genes)
df = t(as.matrix(GetAssayData(so_WT.SGS[feat,], layer = 'data')))


### G1 Phase
DE_wilcox_G1 <- SGS_res_G1$DE %>% arrange(p_val_adj)
feat = DE_wilcox_G1 %>% pull(genes)
df = t(as.matrix(GetAssayData(so_G1[feat,], layer = 'data')))


### S Phase
DE_wilcox_S <- SGS_res_S$DE %>% arrange(p_val_adj)
feat = DE_wilcox_S %>% pull(genes)
df = t(as.matrix(GetAssayData(so_S[feat,], layer = 'data')))

### Saving
write.csv(DE_wilcox_G1,'Results/Kdm6b/Whole/Whole_SGS_G1.csv')
write.csv(DE_wilcox_S,'Results/Kdm6b/Whole/Whole_SGS_S.csv')
write.csv(DE_wilcox_STS_Full,'Results/Kdm6b/Whole/Whole_SGS_Full.csv')

write.csv(SGS_res$HVG_df,'Results/Kdm6b/Whole/Whole_SGS_HVG_Full.csv')

## Invivo
### Pre-processing
so_Whole <- NormalizeData(so_Whole)
res = FindMarkers(so_Whole, ident.1 = 'WT', ident.2 = 'KO', group.by = 'Batch',
                  logfc.threshold = 0.1, min.pct = 0.1, test.use = "wilcox")
DE_wilcox_OG <- res %>% mutate(genes = row.names(res)) %>% arrange(p_val_adj)
head(DE_wilcox_OG)
write.csv(DE_wilcox_OG,'Results/Kdm6b/Whole/Whole_DE_OG.csv')

# Check Seqeuencing depth
so_WT.SGS@meta.data %>% ggplot(aes(x = GoI_SGS,y = nCount_RNA)) +
  geom_boxplot(aes(fill = GoI_SGS), fill = colmap_SGS) +
  ggsignif::geom_signif(comparisons = list(c('Active','Silenced')),
                        map_signif_level = F)+
  theme_classic()

###############################
# Neurons ####
###############################
library(VennDiagram)
library(ggplot2)
library(ggpubr)

SGS_WT = read.csv('Results/Kdm6b/Whole/Whole_SGS_Full.csv')
OG_DE = read.csv('Results/Kdm6b/Whole/Whole_DE_OG.csv')

clipr::write_clip(SGS_WT$genes[1:50])

nSig_SGS = SGS_WT %>% subset(p_val_adj<0.01)
nSig_OG = OG_DE %>% subset(p_val_adj<0.05) %>%
  subset(abs(avg_log2FC)>0.25)

venn1 = venn.diagram(x = list('SGS' = nSig_SGS$genes,
                             'DE' = nSig_OG$genes),
                     filename = NULL, col = 'black', alpha = 0.7,
                     fill = c(hcl.colors(4,'Peach')[2],'white'), cex = 1.5, cat.cex = 1.5,
                     fontfamily = "sans",ext.text = T,ext.dist = -0.1, ext.length = 0.9,
                     cat.fontface = "bold",cat.fontfamily = "sans",
                     cat.dist = c(0.06,0.07))
ggarrange(venn1)
svg('plots/Kdm6b/Venn_Neurons_Kdm6b.svg', height = 4, width = 4)
ggarrange(venn1)+theme(plot.margin = margin(0.5,0.5,0.5,0.5, "in"))
dev.off()

# Random ####
int_rand = c()
for (i in 1:100){
  set.seed(i)
  featn = sample(1:length(rownames(so)), length(SGS_WT$genes[1:200]), replace=F)
  rand_feat = rownames(so)[featn]
  int_rand = append(int_rand,length(intersect(rand_feat,nSig_OG$genes)))
}

int_Mon = length(intersect(SGS_WT$genes[1:200],nSig_OG$genes))
cat_int = data.frame(list(name = c('Random','scSGS'),
                          val = c(mean(int_rand),int_Mon),
                          sd = c(sd(int_rand),sd(int_Mon))))

svg('plots/Kdm6b/Bar_Kdm6b_Random.svg', height = 3, width = 3)
ggplot(cat_int)+
  geom_errorbar( aes(x = name, ymax=val+sd, ymin = val-sd), width=0.3,
                 colour="black", alpha=0.9, linewidth = 0.6)+
  geom_col(aes(x = name,y = val), fill = c('lightgrey', hcl.colors(4,'Peach')[2]),
           width = 0.7, col = 'black') +
  theme_classic(base_size = 15) +
  ylab('Overlap with DE') + xlab(NULL) + NoLegend()
dev.off()

# Pval and Correlation
p2 = SGS_WT %>% filter(p_val_adj < 0.05) %>% filter(genes != 'Kdm6b') %>%
  select(p_val_adj, abs_corr) %>% arrange(p_val_adj) %>%
  mutate(log_pval = -log(p_val_adj)) %>%
  ggplot(aes(x = abs_corr, y = log_pval))+ theme_bw() +
  geom_point(size = 1)+geom_smooth(method = lm, col = 'firebrick2')+
  stat_cor(method="pearson", label.y.npc="top", label.x.npc = "left")+
  ylab('-Log(p-value)')+xlab(expression('abs('*rho*')'))
p2
png('plots/Kdm6b/Correlation_Kdm6b.png', width = 4, height = 4, units = 'in', res = 300)
p2
dev.off()


#############
# Enrichr
#############
library(enrichR)
setEnrichrSite("Enrichr")
dbs <- c("GO_Biological_Process_2023")

## OG DE Whole
gl = nSig_OG$genes[1:200]
OG_Whole_enrichr <- enrichr(gl, dbs)
genenames = sapply(OG_Whole_enrichr[[1]]$Genes[1:100],str_to_title)
genenames = sapply(genenames,
                   function(x) {
                     x = gsub("\\;",", ",as.character(x))
                     gsub("^"," ",as.character(x))
                   })
Ont_OG = OG_Whole_enrichr[[1]] %>% filter(P.value<0.05) %>% pull(Term)

## SGS DE Whole
gl = nSig_SGS$genes[1:200]
SGS_Whole_enrichr <- enrichr(gl, dbs)
genenames = sapply(SGS_Whole_enrichr[[1]]$Genes[1:100],str_to_title)
genenames = sapply(genenames,
                   function(x) {
                     x = gsub("\\;",", ",as.character(x))
                     gsub("^"," ",as.character(x))
                   })

Ont_SGS = SGS_Whole_enrichr[[1]] %>% filter(P.value<0.05) %>% pull(Term)
write.csv(SGS_Whole_enrichr[[1]] %>% filter(P.value<0.05),
          'Results/Kdm6b/Whole/Enrichr_Whole_SGS.csv')

## Pathways Compared
venn = venn.diagram(x = list('SGS' = Ont_SGS,
                             'DE' = Ont_OG),
                    filename = NULL, col = 'black', alpha = 0.7,
                    fill = c(hcl.colors(4,'Peach')[2],'white'), cex = 1.5, cat.cex = 1.5,
                    fontfamily = "sans",ext.text = T,ext.dist = -0.1, ext.length = 0.9,
                    cat.fontface = "bold",cat.fontfamily = "sans",
                    cat.dist = c(0.07,0.07))
svg('plots/Kdm6b/Pathway_Overlaps_Neurons.svg', height = 4, width = 4)
ggarrange(venn) + theme(plot.margin = margin(0.5,0.5,0.5,0.5,'in'))
dev.off()

# Chord diagram
source('SGS/Chord.R')
df = SGS_Whole_enrichr[[1]]
terms = c('Regulation Of Transcription By RNA Polymerase II (GO:0006357)',
          'Regulation Of Neuron Differentiation (GO:0045664)',
          'Regulation Of Double-Strand Break Repair (GO:2000779)')
make_chord(df, terms = terms, font.cex = 1, pal = 'Peach',link.trans =0.3,
           small.gap = 0.2)

# Make the circular plot
svg('plots/Kdm6b/SGS_Kdm6b_Chord.svg', height = 6, width = 6)
make_chord(df, terms = terms, font.cex = 0.8, pal = 'Peach',link.trans = 0.3)
dev.off()


# Unique SGS Genes
Unq_gl = setdiff(nSig_genes_SGS_Full$genes, nSig_genes_OG$genes)[1:200]
Unq_Mon_enrichr <- enrichr(Unq_gl, dbs)
genenames = sapply(Unq_Mon_enrichr[[1]]$Genes[1:50],str_to_title)
genenames = sapply(genenames,
                   function(x) {
                     x = gsub("\\;",", ",as.character(x))
                     gsub("^"," ",as.character(x))
                   })

png('plots/Kdm6b/Neuron_Enrichr_GOBP_Uniq.png', height = 3, width = 11, units = 'in', res = 600 )
plotEnrich(Unq_Mon_enrichr[[1]], showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",
           title = 'GSEA - Unique scSGS responsive genes')+ scale_fill_gradient(low = '#ED90A4',high = '#5AB5E2')+
  geom_text(aes(label = genenames[1:10]), colour = 'white',y = 0, size = 2, hjust = 'left',  fontface = 'bold')
dev.off()

# Common SGS Genes
Comm_gl = intersect(nSig_genes_SGS_Full$genes, nSig_genes_OG$genes)[1:200]
Comm_Mon_enrichr <- enrichr(Comm_gl, dbs)
genenames = sapply(Comm_Mon_enrichr[[1]]$Genes[1:100],str_to_title)
genenames = sapply(genenames,
                   function(x) {
                     x = gsub("\\;",", ",as.character(x))
                     gsub("^"," ",as.character(x))
                   })

png('plots/Kdm6b/Neuron_Enrichr_GOBP_Comm.png', height = 3, width = 10, units = 'in', res = 600 )
plotEnrich(Comm_Mon_enrichr[[1]], showTerms = 10, numChar = 50, y = "Count", orderBy = "Adjusted.P.value",
           title = 'GSEA - Common scSGS responsive genes')+ scale_fill_gradient(low = '#ED90A4',high = '#5AB5E2')+
  geom_text(aes(label = genenames[1:10]), colour = 'white',y = 0, size = 2, hjust = 'left',  fontface = 'bold')
dev.off()
