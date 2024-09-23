futile.logger::flog.threshold(futile.logger::ERROR, name="VennDiagramLogger")
library(Seurat)
library(tidyverse)
library(ggpubr)
library(dittoSeq)
library(VennDiagram)
library(ggplot2)
library(ggsignif)
library(scSGS)

# Set Colormaps
colmap_Batches=c('KO'=hcl.colors(4,'Peach')[2], 'WT'='grey')
colmap_SGS=c('Silenced'=hcl.colors(4,'Peach')[2], 'Active'='grey')

# Load Data ####
so=readRDS('Ccr2/RNA/RNA_glioblastoma_UMAP.rds')
so$Batches=so$sample %>% substr(start=1 ,stop=2)
table(so$Batches)
so$Cell_Overview=recode_factor(so$cluster,
                                 DC1='DC',DC2='DC',DC3='DC',DC4='DC',
                                 'Regulatory T cells'='T cells',
                                 'plasma B cells'='B cells' )
so$Cell_Overview=factor(so$Cell_Overview,
                          levels=c('Mast cells', 'B cells', 'T cells',
                                     'NK cells', 'DC', 'Monocytes',
                                     'prol. TAM', 'TAM 1', 'TAM 2'))

# Cell type bar plot
p1=dittoBarPlot(so, var="Batches", group.by='Cell_Overview',
                  color.panel=colmap_Batches,retain.factor.levels=T,
                  legend.title='Genotype') +ggtitle(NULL) + xlab(NULL)+
  theme_classic(base_size=15)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=0.5))
p1
svg('plots/Ccr2_RNA/Celltype_Batch_Bar.svg', height=4, width=5)
p1
dev.off()

# CTOI
CTOI=so$Cell_Overview
CTOI=c(rep('Other cells',length(so$Cell_Overview)))
CTOI[so$Cell_Overview =='Monocytes']='Monocytes'
so$CTOI=factor(CTOI)

p2=DimPlot(so, reduction="umap", group.by='CTOI', label=T, label.box=T,
             cols=c(hcl.colors(4,'Peach')[2],'lightgrey'),
             label.size=5)+
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        panel.border=element_rect(color='black'))+
  NoLegend()+ggtitle(NULL)
svg('plots/Ccr2_RNA/UMAP_Monocytes.svg', height=4, width=4)
p2
dev.off()

so_mon=subset(so,subset=Cell_Overview =='Monocytes')
so_mon=RunUMAP(so_mon, dims=1:40)
p1=FeaturePlot(subset(so_mon, subset=Batches =='WT'), 'Ccr2', order=T,
                 pt.size=1, cols=c('lightgrey',hcl.colors(4,'Peach')[1]))+
  theme(axis.ticks=element_blank(),
        axis.text=element_blank(),
        axis.line=element_blank(),
        panel.border=element_rect(color='black',linewidth=1))+
  NoLegend() + ggtitle(NULL)
p1
svg('plots/Ccr2_RNA/Ccr2_Monocytes_WT.svg', height=4, width=4)
p1
dev.off()

## Gene Function Analysis ####
Mon_OG_DE=read.csv('Results/Ccr2/Monocytes/Monocytes_SGS_DE_OG.csv', header=T, row.names=1)
Mon_SGS_Full_DE=read.csv('Results/Ccr2/Monocytes/Monocytes_SGS_DE_Full.csv', header=T, row.names=1)
Mon_SGS_G1_DE=read.csv('Results/Ccr2/Monocytes/Monocytes_SGS_DE_G1.csv', header=T, row.names=1)
Mon_SGS_S_DE=read.csv('Results/Ccr2/Monocytes/Monocytes_SGS_DE_S.csv', header=T, row.names=1)

nSig_Mon_OG=Mon_OG_DE %>% subset(p_val_adj<0.05) %>% subset(abs(avg_log2FC) > 0.25)
nSig_Mon_SGS_Full=Mon_SGS_Full_DE %>% subset(p_val_adj<0.01)
nSig_Mon_SGS_G1=Mon_SGS_G1_DE %>% subset(p_val_adj<0.01)
nSig_Mon_SGS_S=Mon_SGS_S_DE %>% subset(p_val_adj<0.01)

# Overlap DE and scSGS
venn2=venn.diagram(x=list(DE=nSig_Mon_OG$genes,
                              SGS=nSig_Mon_SGS_Full$genes[1:200]),
                     filename=NULL, col='black', alpha=0.7,
                     fill=c('white',hcl.colors(4,'Peach')[2]),
                     fontfamily="sans", ext.text=F, cex=1.5, cat.cex=1.5,
                     cat.fontface="bold",cat.fontfamily="sans",
                     cat.dist=c(0.06,0.07))

svg('plots/Ccr2_RNA/Venn_Monocytes_Ccr2.svg', height=4, width=4)
ggarrange(venn2)+theme(plot.margin=margin(0.5,0.5,0.5,0.5, "in"))
dev.off()

# Random genes vs SGS-responsive genes
int_rand=c()
for (i in 1:100){
  set.seed(i)
  featn=sample(1:length(rownames(so)), length(nSig_Mon_SGS_Full$genes[1:200]),
                 replace=F)
  rand_feat=rownames(so)[featn]
  int_rand=append(int_rand,length(intersect(rand_feat,nSig_Mon_OG$genes)))
}

int_Mon=length(intersect(nSig_Mon_SGS_Full$genes[1:200],nSig_Mon_OG$genes))
cat_int=data.frame(list(name=c('Random','scSGS'),
                          val=c(mean(int_rand),int_Mon),
                          sd=c(sd(int_rand),sd(int_Mon))))

svg('plots/Ccr2_RNA/Bar_Monocytes_Ccr2_Random.svg', height=3, width=3)
ggplot(cat_int)+
  geom_errorbar( aes(x=name, ymax=val+sd, ymin=val-sd), width=0.3,
                 colour="black", linewidth=0.6)+
  geom_col(aes(x=name,y=val), width=0.7, col='black',
           fill=c('lightgrey', hcl.colors(4,'Peach')[2])) +
  theme_classic(base_size=15) +
  ylab('Overlap with DE') + xlab(NULL) + NoLegend()
dev.off()


# G1 and S Phase cells vs Full
venn5=venn.diagram(x=list(All=nSig_Mon_SGS_Full$genes,
                              G1=nSig_Mon_SGS_G1$genes,
                              S=nSig_Mon_SGS_S$genes),
                     filename=NULL, col='black',
                     fill=c(hcl.colors(4,'Peach')[2],'white','white'),
                     fontfamily="sans", ext.text=T, cex=1, cat.cex=1.5,
                     cat.fontface="bold",cat.fontfamily="sans",
                     cat.dist=c(0.06,0.09,0.03))

svg('plots/Ccr2_RNA/Venn_Monocytes_Ccr2_Comparison.svg', height=4, width=4)
ggarrange(venn5)+theme(plot.margin=margin(0.5,0.5,0.5,0.5, "in"))
dev.off()

## Enrichr
library(enrichR)
setEnrichrSite("Enrichr")
dbs <- c("GO_Biological_Process_2023")

gl=nSig_Mon_SGS_Full$genes[1:200]
SGS_Mon_enrichr <- enrichr(gl, dbs)
Ont_SGS=SGS_Mon_enrichr[[1]] %>% filter(P.value<0.05) %>% pull(Term)
write.csv(SGS_Mon_enrichr[[1]] %>% filter(P.value<0.05),
          'Results/Ccr2/Monocytes/Enrichr_Monocytes_SGS.csv')

# Compare Pathways
gl=nSig_Mon_OG$genes[1:200]
gl=gl
OG_Mon_enrichr <- enrichr(gl, dbs)
Ont_OG=OG_Mon_enrichr[[1]] %>% filter(P.value<0.05) %>% pull(Term)

venn=venn.diagram(x=list('SGS'=Ont_SGS,
                             'DE'=Ont_OG),
                    filename=NULL, alpha=0.7,
                    fill=c(hcl.colors(4,'Peach')[2],'white'), cex=1.5, cat.cex=1.5,
                    fontfamily="sans",ext.text=F,
                    cat.fontface="bold",cat.fontfamily="sans",
                    cat.dist=c(0.07,0.06))
svg('plots/Ccr2_RNA/Pathway_Overlaps_Monocytes.svg', height=4, width=4)
ggarrange(venn) + theme(plot.margin=margin(0.5,0.5,0.5,0.5,'in'))
dev.off()

## Circos plot - Chord diagram
source('SGS/Chord.R')
df=SGS_Mon_enrichr[[1]]
terms=c('Regulation Of Macrophage Activation (GO:0043030)',
          'Regulation Of Dendritic Cell Differentiation (GO:2001198)',
          'Positive Regulation Of Mononuclear Cell Migration (GO:0071677)',
          'Negative Regulation Of Leukocyte Activation (GO:0002695)',
          'Phagocytosis (GO:0006909)')
make_chord(df, terms=terms, font.cex=1, pal='Peach',link.trans =0.3)

# Make the circular plot
svg('plots/Ccr2_RNA/SGS_Ccr2_Chord.svg', height=6, width=6)
make_chord(df, terms=terms, font.cex=1, pal='Peach',link.trans=0.3)
dev.off()


# Library size Plotting
library(VennDiagram)
library(ggplot2)
library(ggpubr)

SGS_WT=read.csv('Results/Ccr2/Monocytes/Monocytes_SGS_DE_Full.csv', row.names=1) %>%
  arrange(p_val_adj) %>% mutate(Rank=seq.int(nrow(.)))
nSig_SGS=SGS_WT %>% subset(p_val_adj<0.01)

SGS_highlow=read.csv('Results/Ccr2/Monocytes/Monocytes_Highlow_SGS.csv', row.names=1) %>%
  arrange(p_val_adj) %>% mutate(Rank=seq.int(nrow(.)))
nSig_SGS_highlow=SGS_highlow %>% subset(p_val_adj<0.01)

inter=intersect(row.names(SGS_WT),row.names(SGS_highlow))
df_highlow=SGS_WT[inter,] %>% select(p_val_adj,Rank)

venn1=venn.diagram(x=list('SGS'=nSig_SGS$genes[1:200],
                              'Libsize'=nSig_SGS_highlow$genes[1:200]),
                     filename=NULL, col='black', alpha=0.7,
                     fill=c('white',hcl.colors(4,'Peach')[2]), cex=1.5, cat.cex=1.5,
                     fontfamily="sans",ext.text=F,
                     cat.fontface="bold",cat.fontfamily="sans",
                     cat.dist=c(0.06,0.07))
ggarrange(venn1)+theme(plot.margin=margin(0.5,0.5,0.5,0.5, "in"))


