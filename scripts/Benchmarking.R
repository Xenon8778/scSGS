futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
library(Seurat)
library(tidyverse)
library(scSGS)

# Creating Seurat objects and h5ad ####

## Ccr2 Glioblastoma
so = readRDS('Ccr2/RNA/RNA_glioblastoma_UMAP.rds')
so$Batches = so$sample %>% substr(start = 1 ,stop = 2) # Merging Replicates
so_Monocytes = subset(so, subset = cluster == c('Monocytes'))
WT_Mon = subset(so_Monocytes, subset = Batches == "WT")

HVGs = HVG_splinefit(GetAssayData(WT_Mon, layer = 'count'), nHVGs = 5000,
                     dropout.filter = F, diptest = F)
hvgenes = rownames(HVGs %>% filter(HVG == T))


##################
# Benchmarking with STAT1 gene
##################

# 1000 genes 500 cells ####
so_1000_500 = WT_Mon[hvgenes[1:1000],
                     sample(colnames(WT_Mon), size = 500, replace = F)]
saveRDS(so_1000_500, 'Results/Bench/so_1000_500.rds')
cnts = GetAssayData(so_1000_500, layer = 'counts')

# scSGS
gc()
scSGS_out = c()
for (i in 1:5){
  time = system.time(scSGS_res <- scSGS(cnts, GoI = 'Ccr2', calcHVG = T,
                                        nHVG = 500, show.spline = F,
                                        filter_data = F))
  print(time[3])
  scSGS_out = append(scSGS_out,time[3])
}

# scTenifoldKnk
gc()
library(scTenifoldKnk)
Knk_out = c()
for (i in 1:5){
  time = system.time(Knk_res <- scTenifoldKnk(cnts, gKO = 'Ccr2',nc_nNet = 3))
  print(time[3])
  Knk_out = append(Knk_out,time[3])
}

times = data.frame(Knk_out,scSGS_out)

write.csv(scSGS_res$DE,'Results/Bench/SGS_res_1000_500.csv')
write.csv(Knk_res$diffRegulation,'Results/Bench/Knk_res_1000_500.csv')
write.csv(times,'Results/Bench/Knk_SGS_1000_500.csv')

# 1000 genes 1000 cells ####
so_1000_1000 = WT_Mon[hvgenes[1:1000],
                     sample(colnames(WT_Mon), size = 1000, replace = F)]
cnts = GetAssayData(so_1000_1000, layer = 'counts')
saveRDS(so_1000_1000, 'Results/Bench/so_1000_1000.rds')

# scSGS
gc()
scSGS_out = c()
for (i in 1:5){
  time = system.time(scSGS_res <- scSGS(cnts, GoI = 'Ccr2', calcHVG = T,
                                        nHVG = 500, show.spline = F,
                                        filter_data = F))
  print(time[3])
  scSGS_out = append(scSGS_out,time[3])
}

# scTenifoldKnk
gc()
library(scTenifoldKnk)
Knk_out = c()
for (i in 1:5){
  time = system.time(Knk_res <- scTenifoldKnk(cnts, gKO = 'Ccr2',nc_nNet = 3))
  print(time[3])
  Knk_out = append(Knk_out,time[3])
}

times = data.frame(Knk_out,scSGS_out)

write.csv(scSGS_res$DE,'Results/Bench/SGS_res_1000_1000.csv')
write.csv(Knk_res$diffRegulation,'Results/Bench/Knk_res_1000_1000.csv')
write.csv(times,'Results/Bench/Knk_SGS_1000_1000.csv')

# 3000 genes 1000 cells ####
so_3000_1000 = WT_Mon[hvgenes[1:3000],
                      sample(colnames(WT_Mon), size = 1000, replace = F)]
cnts = GetAssayData(so_3000_1000, layer = 'counts')
saveRDS(so_3000_1000, 'Results/Bench/so_3000_1000.rds')


# scSGS
gc()
scSGS_out = c()
for (i in 1:5){
  time = system.time(scSGS_res <- scSGS(cnts, GoI = 'Ccr2', calcHVG = T,
                                        nHVG = 500, show.spline = F,
                                        filter_data = F))
  print(time[3])
  scSGS_out = append(scSGS_out,time[3])
}

# scTenifoldKnk
gc()
library(scTenifoldKnk)
Knk_out = c()
for (i in 1:5){
  time = system.time(Knk_res <- scTenifoldKnk(cnts, gKO = 'Ccr2',nc_nNet = 3))
  print(time[3])
  Knk_out = append(Knk_out,time[3])
}

times = data.frame(Knk_out,scSGS_out)

write.csv(scSGS_res$DE,'Results/Bench/SGS_res_3000_1000.csv')
write.csv(Knk_res$diffRegulation,'Results/Bench/Knk_res_3000_1000.csv')
write.csv(times,'Results/Bench/Knk_SGS_3000_1000.csv')

# 5000 genes 1000 cells ####
so_5000_1000 = WT_Mon[hvgenes,
                      sample(colnames(WT_Mon), size = 1000, replace = F)]
cnts = GetAssayData(so_5000_1000, layer = 'counts')
saveRDS(so_5000_1000, 'Results/Bench/so_5000_1000.rds')

# scSGS
gc()
scSGS_out = c()
for (i in 1:5){
  time = system.time(scSGS_res <- scSGS(cnts, GoI = 'Ccr2', calcHVG = T,
                                        nHVG = 500, show.spline = F,
                                        filter_data = F))
  print(time[3])
  scSGS_out = append(scSGS_out,time[3])
}

# scTenifoldKnk
gc()
library(scTenifoldKnk)
Knk_out = c()
for (i in 1:5){
  time = system.time(Knk_res <- scTenifoldKnk(cnts, gKO = 'Ccr2',nc_nNet = 3))
  print(time[3])
  Knk_out = append(Knk_out,time[3])
}

times = data.frame(Knk_out,scSGS_out)

write.csv(scSGS_res$DE,'Results/Bench/SGS_res_5000_1000.csv')
write.csv(Knk_res$diffRegulation,'Results/Bench/Knk_res_5000_1000.csv')
write.csv(times,'Results/Bench/Knk_SGS_5000_1000.csv')

# 5000 genes 3000 cells ####
so_5000_3000 = WT_Mon[hvgenes,
                      sample(colnames(WT_Mon), size = 3000, replace = F)]
cnts = GetAssayData(so_5000_3000, layer = 'counts')
saveRDS(so_5000_3000, 'Results/Bench/so_5000_3000.rds')

# scSGS
gc()
scSGS_out = c()
for (i in 1:5){
  time = system.time(scSGS_res <- scSGS(cnts, GoI = 'Ccr2', calcHVG = T,
                                        nHVG = 500, show.spline = F,
                                        filter_data = F))
  print(time[3])
  scSGS_out = append(scSGS_out,time[3])
}

# scTenifoldKnk
gc()
library(scTenifoldKnk)
Knk_out = c()
for (i in 1:5){
  time = system.time(Knk_res <- scTenifoldKnk(cnts, gKO = 'Ccr2',nc_nNet = 3))
  print(time[3])
  Knk_out = append(Knk_out,time[3])
}

times = data.frame(Knk_out,scSGS_out)

write.csv(scSGS_res$DE,'Results/Bench/SGS_res_5000_3000.csv')
write.csv(Knk_res$diffRegulation,'Results/Bench/Knk_res_5000_3000.csv')
write.csv(times,'Results/Bench/Knk_SGS_5000_3000.csv')

# Diffrential Expression ####
set.seed(0)
so_OG = so_Monocytes[hvgenes,
                     sample(colnames(so_Monocytes), size = 3000, replace = F)]
saveRDS(so_OG,'Results/Bench/so_OG.rds')

DE_res = FindMarkers(so_OG, ident.1 = 'KO', ident.2 = 'WT', group.by = 'Batches')
write.csv(DE_res,'Results/Bench/DE_res_5000_3000.csv')

### Plotting Times ####
library(ggplot2)
library(dplyr)
library(ggpubr)

times = read.csv('Results/Bench/Knk_SGS_1000_500.csv', row.names = 1)

df_1000_500 <- data.frame(
  name=colnames(times),
  value=colMeans(times),
  sd=colSds(as.matrix(times))
) %>% mutate(dataset = '1000_500') %>% mutate(name = c('Knk','SGS'))
df_1000_500

times = read.csv('Results/Bench/Knk_SGS_1000_1000.csv', row.names = 1)
df_1000_1000 <- data.frame(
  name=colnames(times),
  value=colMeans(times),
  sd=colSds(as.matrix(times))
) %>% mutate(dataset = '1000_1000') %>% mutate(name = c('Knk','SGS'))
df_1000_1000

times = read.csv('Results/Bench/Knk_SGS_3000_1000.csv', row.names = 1)
df_3000_1000 <- data.frame(
  name=colnames(times),
  value=colMeans(times),
  sd=colSds(as.matrix(times))
) %>% mutate(dataset = '3000_1000') %>% mutate(name = c('Knk','SGS'))
df_3000_1000

times = read.csv('Results/Bench/Knk_SGS_5000_1000.csv', row.names = 1)
df_5000_1000 <- data.frame(
  name=colnames(times),
  value=colMeans(times),
  sd=colSds(as.matrix(times))
) %>% mutate(dataset = '5000_1000') %>% mutate(name = c('Knk','SGS'))
df_5000_1000

times = read.csv('Results/Bench/Knk_SGS_5000_3000.csv', row.names = 1)
df_5000_3000 <- data.frame(
  name=colnames(times),
  value=colMeans(times),
  sd=colSds(as.matrix(times))
) %>% mutate(dataset = '5000_3000') %>% mutate(name = c('Knk','SGS'))
df_5000_3000

df_times = rbind(df_1000_500,df_1000_1000,df_3000_1000,df_5000_1000,df_5000_3000)
df_times$name = factor(df_times$name, levels = c('Knk','SGS'))
write.csv(df_times,'Results/Bench/All_times.csv')

df_times = read.csv('Results/Bench/All_times.csv')

png('plots/Bench Runtimes.png', width = 5, height = 4, res = 300, units = 'in')
ggplot(df_times, aes(x=dataset, y=value, group = name)) +
  geom_errorbar(aes(ymin=value-sd, ymax=value+sd),
                width=0.4, colour="black", linewidth=0.5,
                position = position_dodge(width = 0.9)) +
  geom_bar(stat="identity", aes(fill = name),
           position = position_dodge(width = 0.9))+
  scale_fill_manual(values = hcl.colors(3,'Set2')) +
  theme_classic() + xlab('Dataset (Genes_Cells)') + ylab('Runtime (seconds)') +
  labs(fill = 'Tool')
dev.off()

### Plotting Overlaps ####
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
library(VennDiagram)

DE_res = read.csv('Results/Bench/DE_res_5000_3000.csv', row.names = 1) %>%
  mutate(genes = rownames(.))
SGS_res = read.csv('Results/Bench/SGS_res_5000_3000.csv')
Knk_res = read.csv('Results/Bench/Knk_res_5000_3000.csv')
GenKI_res = read.csv('Results/Bench/GenKI_res_5000_3000.csv')

sigDE = DE_res %>% arrange(p_val_adj) %>% filter(p_val_adj < 0.05) %>% pull(genes)
sigSGS = SGS_res %>% arrange(p_val_adj) %>% head(100) %>% pull(genes)
sigKnk = Knk_res %>% arrange(p.value) %>% head(100) %>% pull(gene)
sigGenKI = GenKI_res %>% arrange(rank) %>% head(100) %>% pull(X)

venn1 = venn.diagram(list(SGS = sigSGS,
                          DE = sigDE,
                          scTenifoldknk = sigKnk,
                          GenKI = sigGenKI),
                     filename = NULL,
                     fill = c(hcl.colors(4,'Peach')[2],'white','white','white'),
                     cat.fontfamily = 'sans', fontfamily = 'sans', margin = 0.05,
                     cat.fontface = "bold",
                     cat.dist = c(0.07,0.07,0.1,0.09))
ggarrange(venn1)

svg('plots/Bench Overlap.svg', width = 4, height = 4)
ggarrange(venn1)
dev.off()
