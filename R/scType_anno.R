## Load Libs
library(Seurat)
library(openxlsx)
library(Matrix)
library(dplyr)
library(stringr)

prep_panglaodb <- function(tissue = 'All', species = 'Mm'){

  # Download and read data
  url <- "https://raw.githubusercontent.com/Xenon8778/Auto_cell_annot/main/PanglaoDB_markers_27_Mar_2020.tsv"
  db <- read.delim(url, sep = "\t")

  # Subset relevant columns
  db_sub <- db %>%
    select(species, `official.gene.symbol`, `cell.type`, `ubiquitousness.index`,`organ`)

  # Select species
  species <- species # Replace with your desired species

  if ('All' %in% tissue){
    db_sub <- db_sub %>%
      filter(species %in% c(species, "Mm Hs", "Hs Mm")) %>%
      filter(`ubiquitousness.index` < 0.1)
  } else {
    db_sub <- db_sub %>%
      filter(species %in% c(species, "Mm Hs", "Hs Mm")) %>%
      filter(organ == tissue) %>%
      filter(`ubiquitousness.index` < 0.1)
  }

  # Prepare data by cell type
  db_prep <- db_sub %>%
    group_by(cell.type,organ) %>%
    reframe(geneSymbolmore = list(`official.gene.symbol`))
  db_prep$geneSymbolmore1 = sapply(db_prep$geneSymbolmore, paste, collapse=",")
  db_prep <- db_prep %>%
    mutate(tissueType = organ) %>% mutate(cellName = cell.type) %>%
    mutate(geneSymbolmore2 = NA) %>% mutate(shortName = cell.type) %>% ungroup() %>%
    select(tissueType,cellName,geneSymbolmore1,geneSymbolmore2,shortName)

  return(db_prep)
}

gene_sets_prepare <- function(db_file, cell_type){

  cell_markers = db_file
  if ("All" %in% cell_type){
    cell_markers = cell_markers
  } else {
    cell_markers = cell_markers[cell_markers$tissueType == cell_type,]
  }
  cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1); cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)

  # correct gene symbols from the given DB (up-genes)
  cell_markers$geneSymbolmore1 = sapply(1:nrow(cell_markers), function(i){

    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore1[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)

    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })

  # correct gene symbols from the given DB (down-genes)
  cell_markers$geneSymbolmore2 = sapply(1:nrow(cell_markers), function(i){

    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$geneSymbolmore2[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)

    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(checkGeneSymbols(markers_all)$Suggested.Symbol))})
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })

  cell_markers$geneSymbolmore1 = gsub("///",",",cell_markers$geneSymbolmore1);cell_markers$geneSymbolmore1 = gsub(" ","",cell_markers$geneSymbolmore1)
  cell_markers$geneSymbolmore2 = gsub("///",",",cell_markers$geneSymbolmore2);cell_markers$geneSymbolmore2 = gsub(" ","",cell_markers$geneSymbolmore2)

  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore1[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$geneSymbolmore2[j]),",")))); names(gs2) = cell_markers$cellName

  list(gs_positive = gs, gs_negative = gs2)
}

sctype_score <- function(scRNAseqData, scaled = !0, gs, gs2 = NULL, gene_names_to_uppercase = !0, ...){

  # check input matrix
  if(!is.matrix(scRNAseqData)){
    warning("scRNAseqData doesn't seem to be a matrix")
  } else {
    if(sum(dim(scRNAseqData))==0){
      warning("The dimension of input scRNAseqData matrix equals to 0, is it an empty matrix?")
    }
  }

  # marker sensitivity
  marker_stat = sort(table(unlist(gs)), decreasing = T);
  marker_sensitivity = data.frame(score_marker_sensitivity = scales::rescale(as.numeric(marker_stat), to = c(0,1), from = c(length(gs),1)),
                                  gene_ = names(marker_stat), stringsAsFactors = !1)

  # convert gene names to Uppercase
  if(gene_names_to_uppercase){
    rownames(scRNAseqData) = toupper(rownames(scRNAseqData));
  }

  # subselect genes only found in data
  names_gs_cp = names(gs); names_gs_2_cp = names(gs2);
  gs = lapply(1:length(gs), function(d_){
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  gs2 = lapply(1:length(gs2), function(d_){
    GeneIndToKeep = rownames(scRNAseqData) %in% as.character(gs2[[d_]]); rownames(scRNAseqData)[GeneIndToKeep]})
  names(gs) = names_gs_cp; names(gs2) = names_gs_2_cp;
  cell_markers_genes_score = marker_sensitivity[marker_sensitivity$gene_ %in% unique(unlist(gs)),]

  # z-scale if not
  if(!scaled) Z <- t(scale(t(scRNAseqData))) else Z <- scRNAseqData

  pb = txtProgressBar(min = 0,
                      max = nrow(cell_markers_genes_score),
                      initial = 0,
                      style = 3)

  # multiple by marker sensitivity
  for(jj in 1:nrow(cell_markers_genes_score)){
    Z[cell_markers_genes_score[jj,"gene_"], ] = Z[cell_markers_genes_score[jj,"gene_"], ] * cell_markers_genes_score[jj, "score_marker_sensitivity"]
    setTxtProgressBar(pb,jj)
  }
  close(pb)

  # subselect only with marker genes
  Z = Z[unique(c(unlist(gs),unlist(gs2))), ]

  # combine scores
  es = do.call("rbind", lapply(names(gs), function(gss_){
    sapply(1:ncol(Z), function(j) {
      gs_z = Z[gs[[gss_]], j]; gz_2 = Z[gs2[[gss_]], j] * -1
      sum_t1 = (sum(gs_z) / sqrt(length(gs_z))); sum_t2 = sum(gz_2) / sqrt(length(gz_2));
      if(is.na(sum_t2)){
        sum_t2 = 0;
      }
      sum_t1 + sum_t2
    })
  }))

  dimnames(es) = list(names(gs), colnames(Z))
  es.max <- es[!apply(is.na(es) | es == "", 1, all),] # remove na rows

  es.max
}

Annotate_cells <- function(data = x, res = 0.8, tissue = "All", Scale = F, annot_only = F,
                           filter_conf = T, ndims = 50, DB = 'scType' , species = 'Mm'){
  #> data - is a count matrix
  #> res - resolution for graph based clustering
  #> tissue - Tissue type for which we are looking markers for
  #> DB - either 'scType' or 'PanDB'

  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  if (DB == 'PanDB'){
    print("Using PanglaoDB")
    dataDB = prep_panglaodb(tissue = tissue, species = species)
  } else {
    print("Using ScType")
    dataDB = read.xlsx(db_)
    if ("All" %in% tissue){
      print("Checking across all celltypes.\n Specify following cell types for in-depth analysis")
      print(unique(dataDB$tissueType))
    }
  }

  data = CreateSeuratObject(CreateAssayObject(data))
  data <- NormalizeData(data)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  if (isTRUE(Scale)){data <- ScaleData(data, features = rownames(data))}
  data <- RunPCA(data,verbose = F)
  data <- FindNeighbors(data, dims = 1:ndims)
  data <- FindClusters(data, resolution = res)
  data <- RunUMAP(data, dims = 1:ndims)

  # load libraries and functions
  lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)


  # prepare gene sets
  gs_list = gene_sets_prepare(dataDB,tissue)
  print('Genesets Prepared')

  # get cell-type by cell matrix
  es.max = sctype_score(scRNAseqData = data[["RNA"]]@scale.data, scaled = TRUE,
                        gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  print('Scoring Done')

  # merge by cluster
  cL_resutls = do.call("rbind", lapply(unique(data@meta.data$seurat_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(data@meta.data[data@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(data@meta.data$seurat_clusters==cl)), 10)
  }))
  sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)

  # set low-confident (low ScType score) clusters to "unknown"
  if (filter_conf == T){
    sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
    }
  print(head(sctype_scores[,1:3]))

  data@meta.data$scType_anno = ""
  for(j in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==j,];
    data@meta.data$scType_anno[data@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
  }
  if (annot_only == T){
    return(unlist(data$scType_anno))
  }
  else{return(data)}
}

#> library(ggplot2)
#> library(dittoSeq)
#> library(gridExtra)
#> DimPlot(so, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'scType_anno')
#> dittoBarPlot(so,var = "Batches", group.by = "scType_anno")+
#>  scale_fill_manual("Batches",values=c("#F68282", "#7CAE00", "#00BFC4","#C77CFF"))
