## Load Libs
library(circlize)
library(dplyr)

make_chord = function(df, terms, pal = 'Set2', link.wt = 1, species = 'Mm',
                      link.trans = 0.3, font.cex = 0.8,
                      small.gap = 1, big.gap = 5){
  # Load the circlize library
  library(circlize)

  df1 = df %>% filter(Term %in% terms)
  row.names(df1) = df1$Term
  df1 = df1[terms,]
  df1

  genes = str_split(unlist(df1['Genes']),';')
  allg = unique(unlist(genes))


  short = lapply(str_split(terms,' '),function(x) tail(x,1))
  short = str_replace(short,':','')
  short = str_replace(short,'[)]','')
  short = str_replace(short,'[(]','')
  short

  data = list()
  for (i in allg){
    l = c()
    for (j in 1:length(terms)){
      print(j)
      if(i %in% unlist(genes[j])){
        l = append(l,0.8)
      }else{
        l = append(l,0)
      }
    }
    data[i] = list(l)
  }

  print(short)
  data
  data = as.data.frame(data, row.names = short)
  if (species == 'Mm'){
    colnames(data) = str_to_title(colnames(data))
  }
  data = data %>% mutate(terms = row.names(.))
  data = melt(data, id = 'terms')

  cols = c(hcl.colors(length(terms),pal),rep(hcl.colors(1,pal),length(allg)))
  names(cols) = c(short,as.character(unique(data$variable)))

  circos.par(start.degree = 180)
  chordDiagram(data, annotationTrack = "grid", preAllocateTracks = 2,
               grid.col = cols, small.gap = small.gap, big.gap = big.gap,
               transparency = link.trans,
               col = hcl.colors(length(terms),pal),
               link.border = alpha(1,0.3), link.lwd = link.wt)
  circos.trackPlotRegion(track.index = 2, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + 0, sector.name,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                cex = font.cex)
  }, bg.border = NA)
  circos.clear()
}

# Use case
# make_chord(df, terms, pal = 'Peach')
