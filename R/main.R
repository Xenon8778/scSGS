#' @title Spline-HVG
#' @description Compute Highly Variable Genes
#' @export HVG_splinefit
#' @import dplyr
#' @import plotly
#' @rawNamespace import(ggplot2, except = last_plot)
#' @importFrom sparseMatrixStats rowSds
#' @importFrom Matrix rowSums rowMeans
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Biobase matchpt
#' @importFrom stats p.adjust predict quantile smooth.spline
#' @importFrom methods as
#' @author Shreyan Gupta <xenon8778@tamu.edu>
#' Spline HVG function
#' This function computed highly variable genes from an scRNAseq count matrix.

HVG_splinefit <- function(adata, degf = 15,
                          spar = 0.75, nHVGs = 2000, show.spline = FALSE,
                          use.ndist = T){

  ## Gene Statistics
  print('Computing Gene Statistics')
  bin_data <- adata == 0
  Dropout = Matrix::rowSums(bin_data)
  Dropout = Dropout/ncol(adata)
  Means = Matrix::rowMeans(adata)
  SDs = sparseMatrixStats::rowSds(adata)
  CV = SDs/Means
  splinefit_df = as.data.frame(Means)
  splinefit_df$CV = CV
  splinefit_df$Dropout = Dropout
  splinefit_df$logMean = log(Means+1)
  splinefit_df$logCV = log(CV+1)
  splinefit_df = splinefit_df %>% arrange(Dropout)
  splinefit_df = splinefit_df %>% filter(Dropout>0.01 & Dropout<0.99) # Trimming ends

  # Calculate the differences and the squared differences between consecutive elements
  diff_lgu <- diff(splinefit_df$logMean)
  diff_lgcv <- diff(splinefit_df$logCV)
  diff_dropr <- diff(splinefit_df$Dropout)
  diff_squared_sum <- diff_lgu^2 + diff_lgcv^2 + diff_dropr^2

  # Calculate the cumulative sum
  s <- c(0, cumsum(sqrt(diff_squared_sum)))
  xyz = splinefit_df %>% select(logMean,logCV,Dropout)

  fitx <- smooth.spline(s, splinefit_df$logMean, df = degf, spar = spar)
  fity <- smooth.spline(s, splinefit_df$logCV, df = degf, spar = spar)
  fitz <- smooth.spline(s, splinefit_df$Dropout, df = degf, spar = spar)

  # Computing Distances
  print('Computing Distances')
  xyz1 = cbind(predict(fitx,s)$y,predict(fity,s)$y,predict(fitz,s)$y)
  xyz1 = as.data.frame(xyz1)
  colnames(xyz1) = c('logMean','logCV','Dropout')
  rownames(xyz1) = rownames(splinefit_df)
  euclidean <- function(a, b) sqrt(sum((a - b)^2))
  if (use.ndist == F){
    Dist_HVG = c()
    pb = txtProgressBar(min = 0, max = nrow(xyz), initial = 0, style = 3)
    for (i in 1:nrow(xyz)){
      Dist_HVG = append(Dist_HVG,euclidean(xyz[i,],xyz1[i,]))
      setTxtProgressBar(pb,i)
    }
    close(pb)
  } else {
    Dist_HVG = c()
    df = as.matrix(xyz1)
    pb = txtProgressBar(min = 0, max = nrow(xyz), initial = 0, style = 3)
    for (i in 1:nrow(xyz)){
      p1 = as.matrix(xyz[i,])
      near_point = matchpt(p1, df)
      Dist_HVG = append(Dist_HVG,near_point$distance)
      setTxtProgressBar(pb,i)
    }
    close(pb)
  }
  splinefit_df$Distance = Dist_HVG
  splinefit_df$splinex = xyz1$logMean
  splinefit_df$spliney = xyz1$logCV
  splinefit_df$splinez = xyz1$Dropout
  top_HVG = splinefit_df %>% top_n(n = nHVGs, wt = Distance)
  mask = rownames(splinefit_df) %in% rownames(top_HVG)
  splinefit_df$HVG = mask
  splinefit_df = splinefit_df %>% arrange(Dropout)

  # Plotting in 3D
  if(show.spline == TRUE){
    fig = plot_ly(splinefit_df, x = ~logMean, y = ~logCV, z = ~Dropout, mode = 'lines+marker',
                  color = ~HVG, colors = c('grey','red'))
    fig <- fig %>% add_markers(mode = 'marker', type = "scatter3d", opacity = 0.5,
                               marker = list(size = 2))
    fig <- fig %>% add_trace(splinefit_df, x = ~splinex, y = ~spliney, z = ~splinez, type="scatter3d", mode = 'lines+marker',
                             opacity = 1, line = list(width = 5, color = 'black', reverscale = FALSE))
    fig <- fig %>% layout(scene = list(xaxis = list(title = 'log(Mean)'),
                                       yaxis = list(title = 'log(CV)'),
                                       zaxis = list(title = 'Dropout rate')))
    print(fig)
  }


  splinefit_df = splinefit_df[order(splinefit_df$Distance, decreasing = T),]
  return(splinefit_df)
}

#' @export scSGS
#' @importFrom sparseMatrixStats rowSums2 colSums2
#' @import dplyr
#' @import openxlsx
#' @import Seurat
#' @import presto
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Matrix rowSums rowMeans colSums
#' @importFrom methods is
#' @title scSGS
#' @author Shreyan Gupta <xenon8778@tamu.edu>
#' @description Predict SGS-reponsive genes
#' @return A result list with - DV Genes dataframe, SGS-DE dataframe, TF List

scSGS <- function(data, GoI, nHVG = 500, HVG_algo = 'splinefit',
                  show.spline = FALSE, calcHVG = FALSE,
                  nfeatures = 200, ncells = 10, norm_sep = F,
                  rm.mt = FALSE, rm.rp = FALSE, filter_data = TRUE){

  #> data -> Gene x Cell count matrix in the form of a data frame, matrix or sparse matrix with
  #>         rows as genes and columns as cells. Row names should be gene names.
  #> GoI -> Gene of Interest for scSGS analysis

  # data = mat
  # GoI = 'STAT1'
  # HVG_algo = 'splinefit'
  # nHVG = 500
  # show.spline = FALSE

  # Converting to Sparse matrix
  if (is(class(data), "dgCMatrix")){
    print('Converting to sparse')
    spmat = as.matrix(data)
    spmat = as(spmat, "dgCMatrix")
  } else {
    spmat = data
  }

  res  = list()

  # Making binary counts
  GoI = GoI
  GoI_Mask = c()
  for (i in spmat[GoI,]){
    if (i > 0){GoI_Mask=append(GoI_Mask,"Active")
    } else {GoI_Mask=append(GoI_Mask,"Silenced")}
  }
  print("SGS counts:")
  print(table(GoI_Mask))
  GoI_Mask = as.factor(GoI_Mask)
  res$GoI_Mask = GoI_Mask

  # Filtering
  if (filter_data == TRUE){
    Cell_Mask <- sparseMatrixStats::colSums2(spmat > 0) > nfeatures &
      sparseMatrixStats::colSums2(spmat > 0) < quantile(sparseMatrixStats::colSums2(spmat > 0),0.99)
    spmat <- spmat[,Cell_Mask]
    spmat <- spmat[sparseMatrixStats::rowSums2(spmat) > ncells,]
  }

  # Making binary counts of filtered data
  GoI = GoI
  GoI_Mask = c()
  for (i in spmat[GoI,]){
    if (i > 0){GoI_Mask=append(GoI_Mask,"Active")
    } else {GoI_Mask=append(GoI_Mask,"Silenced")}
  }
  print("Filtered counts:")
  print(table(GoI_Mask))
  GoI_Mask = as.factor(GoI_Mask)

  # Filter Mitochondrial and Ribosomal protein genes
  if (rm.rp == T){
    spmat <- spmat[!startsWith(rownames(spmat),'Rps'), ]
    spmat <- spmat[!startsWith(rownames(spmat),'Rpl'), ]
    spmat <- spmat[!startsWith(rownames(spmat),'RPS'), ]
    spmat <- spmat[!startsWith(rownames(spmat),'RPL'), ]
  }
  if (rm.mt == T){
    spmat <- spmat[!startsWith(rownames(spmat),'mt-'), ]
    spmat <- spmat[!startsWith(rownames(spmat),'MT-'), ]
  }

  # Normalizing Data Separately - Log-Normalization and scaling to 10,000
  if (norm_sep == F){
    matnorm = spmat
    matnorm@x = matnorm@x / rep.int(Matrix::colSums(matnorm), diff(matnorm@p))
    matnorm@x = log(matnorm@x*10000+1) #Scale factor 10000
  } else {
    matnorm = spmat
  }

  # Compute Highly Variable Genes
  if (calcHVG == T){
    print("Computing Highly Variable Genes")
    # Identifying HVGs with Spline-HVG
    if (HVG_algo == 'splinefit'){
      HVG_df = HVG_splinefit(spmat, nHVGs = nHVG, show.spline = show.spline)
      res$HVG_df = HVG_df
      res$HV_genes = row.names(HVG_df[HVG_df$HVG == T &
                                        HVG_df$Dropout>0.25 &
                                        HVG_df$Dropout<0.75,])
    }

    # Identifying HVGs with Seurat
    if (HVG_algo == 'seurat'){
      seuratmat = CreateAssayObject(matnorm)
      seuratmat = FindVariableFeatures(seuratmat, selection.method = "vst", nfeatures = 500)
      HVG_df = seuratmat@meta.features
      res$HVG_df = HVG_df
      res$HV_genes = row.names(HVG_df[HVG_df$vst.variable == T,])
    }

    ## Load TF Database to check if gene is a known TF
    Tf_db = read.xlsx('https://static-content.springer.com/esm/art%3A10.1038%2Fs41587-020-0742-6/MediaObjects/41587_2020_742_MOESM3_ESM.xlsx')[,1]
    res$TF = res$HV_genes[res$HV_genes %in% Tf_db]

    ## Check if GoI is HVG
    if (GoI %in% res$HV_genes == T){
      print("The Gene of Interest is an HVG")
    } else {
      print("Beware! The Gene of Interest is not variable enough.")
    }
  }

  # Splitting the matrices
  SGS_WT = matnorm[,GoI_Mask=='Active']
  SGS_KO = matnorm[,GoI_Mask=='Silenced']

  if (norm_sep == T){
    # Normalizing SGS_WT - Log-Normalization and scaling to 10,000
    SGS_WT@x = SGS_WT@x / rep.int(colSums(SGS_WT), diff(SGS_WT@p))
    SGS_WT@x = log(SGS_WT@x*10000+1) #Scale factor 10000

    # Normalizing SGS_KO - Log-Normalization and scaling to 10,000
    SGS_KO@x = SGS_KO@x / rep.int(colSums(SGS_KO), diff(SGS_KO@p))
    SGS_KO@x = log(SGS_KO@x*10000+1) #Scale factor 10000
  }

  # Calculate percent expressed
  thresh.min <- 0
  pct.KO <- round(
    x = rowSums(SGS_KO > thresh.min) /
      dim(SGS_KO)[2],
    digits = 3
  )
  pct.WT <- round(
    x = rowSums(x = SGS_WT > thresh.min) /
      dim(SGS_WT)[2],
    digits = 3
  )

  # Calculate fold change
  FC_KO <- log(x = (rowSums(SGS_KO) + 1)/ncol(SGS_WT), base = 2)
  FC_WT <- log(x = (rowSums(SGS_WT) + 1)/ncol(SGS_WT), base = 2)
  fc <- (FC_KO - FC_WT)
  fc.results <- as.data.frame(x = cbind(fc, pct.KO, pct.WT))
  colnames(fc.results) <- c('avg_log2FC', "pct.KO", "pct.WT")

  # Filtering based on pct
  alpha.min = pmax(fc.results$pct.KO, fc.results$pct.WT)
  names(x = alpha.min) = rownames(x = fc.results)
  features = names(x = which(x = alpha.min >= 0.1))
  fc.results = fc.results[features,]

  # Filtering based on FC
  beta.min <- abs(fc.results$avg_log2FC)
  names(x = beta.min) <- rownames(x = fc.results)
  features <- names(x = which(x = beta.min >= 0.1))
  fc.results = fc.results[features,]

  # Wilcoxon Rank-sum Test
  DE = wilcoxauc(matnorm[features,], GoI_Mask)
  DE = DE %>% filter(group == 'Silenced')
  DE.res = fc.results
  DE.res$p_val = DE$pval
  DE.res$p_val_adj = p.adjust(p = DE.res$p_val,
                              method = "BH",
                              n = nrow(matnorm[features,]))
  DE.res = DE.res %>% arrange(p_val_adj) %>% mutate(genes = row.names(.))

  res$DE = DE.res
  return(res)
}

get.DE <- function(res){
  return(res$DE) # Returns the Differential Expression dataframe
}

get.HVGenes <- function(res){
  return(res$HV_genes) #Returns the n Highly variable genes
}


#' @export get.Enrichment
#' @import enrichR
#' @title get.Enrichment
#' @author Shreyan Gupta <xenon8778@tamu.edu>
#' @description Get enrichR gene function results
#' @return A result list with - enrichR DF, anrichr barplot

get.Enrichment <- function(res, db = "GO_Biological_Process_2023", ngenes = NULL,
                           plot.it = T){

  #> Performs gene set enrichment analysis with Enrichr
  #>
  #> Chen EY, Tan CM, Kou Y, Duan Q, Wang Z, Meirelles GV, Clark NR, Ma'ayan A.
  #> Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool.
  #> BMC Bioinformatics. 2013;128(14)
  #>
  #> Kuleshov MV, Jones MR, Rouillard AD, Fernandez NF, Duan Q, Wang Z, Koplev S,
  #> Jenkins SL, Jagodnik KM, Lachmann A, McDermott MG, Monteiro CD, Gundersen GW, Ma'ayan A.
  #> Enrichr: a comprehensive gene set enrichment analysis web server 2016 update.
  #> Nucleic Acids Research. 2016; gkw377.
  #>
  #> res -> scSGS result
  #> db -> Pathway database to use
  #> ngenes -> number of significant genes to use, else uses all significant genes

  # Load Enrichr
  setEnrichrSite("Enrichr")
  dbs <- c(db)

  # Extracting Significant Gene list
  if (is.null(ngenes)){
    DE_df = res$DE %>% arrange(p_val) %>% filter(p_val < 0.05)
    DE_gl = DE_df %>% pull(genes)
  } else {
    DE_df = res$DE %>% arrange(p_val) %>% filter(p_val < 0.05)
    DE_gl = DE_df[1:ngenes,] %>% pull(genes)
  }

  # Performing Enrichment with Enrichr
  enrichr_df <- enrichr(DE_gl, db)

  # plotting
  if (plot.it == T){
    p1 = plotEnrich(enrichr_df[[1]], showTerms = 10, numChar = 100,
                    y = "Count", orderBy = "Adjusted.P.value", title = db)+
      scale_fill_gradient(low = '#ED90A4',high = '#5AB5E2')
    plot(p1)
  }
  return(enrichr_df)
}


#' # Test Code
#' pbmc10k = readRDS('PBMC/pbmc10k_annot.rds')
#' mat = GetAssayData(subset(pbmc10k, subset = scType_anno == 'Naive CD4+ T cells'),
#'                    layer = 'counts')
#' mat = as.data.frame(mat)
#' res = scSGS(mat, GoI = 'STAT1', nHVG = 500, rm.mt = T, rm.rp = T, HVG_algo = 'splinefit')
#' get.DE(res)
#' get.Enrichment(res, db = "GO_Biological_Process_2023", ngenes = 200)
#' get.STRINGdb(res, version = "12.0", species = 9606, ngenes = 30)
