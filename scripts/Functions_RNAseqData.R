##########################################################################
##########################################################################
# Project: NT organoid patterning 
# Script purpose: a collection of functions used for Sequencing data analysis 
# - Elena's RNA-seq time series
# - Foxa2 positive and negative at Day5
# - TM3 data at Day3
# - Human organoid single-cell RNA-seq data from Adran
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Jul 16 09:42:08 2021
##########################################################################
##########################################################################
##########################################
# edit a list of genes in FGF, BMP and Wnt 
##########################################
Manual.curate.geneList.signalingPathways = function()
{
  # examples = unique(c('Foxa2', # FoxA 
  #                     'Lef1', 'Mapk1', rownames(fpm)[grep('Smad', rownames(fpm))], 
  #                     rownames(fpm)[grep('Wnt|Dkk|Tcf', rownames(fpm))],
  #                     rownames(fpm)[grep('Bmp', rownames(fpm))], 'Nog', 'Chrd', 'Runx1', 'Runx2',  'Smad6', 'Id1', 'Id3',
  #                     rownames(fpm)[grep('Acvr', rownames(fpm))],
  #                     rownames(fpm)[grep('Fgf', rownames(fpm))], 
  #                     'Dusp1', 'Dusp10', 'Dusp27', 'Dusp4', 'Dusp5', 'Mapk10', 'Mapk4', 'Mapk8ip2', 'Spry4', 'Rbpj', 'Hes1', 'Hes5',
  #                     'Hes7', 'Hey1', 'Hey2',
  #                     rownames(fpm)[grep('Notch|Jag|Dll|Dlk', rownames(fpm))], 
  #                     rownames(fpm)[grep('Tgf', rownames(fpm))]
  # ))
  
  require(openxlsx)
  xx.files = list.files(path = '../data', pattern =  'xlsx', full.names = TRUE)
  xx.files = xx.files[grep('GO_term', xx.files, invert = FALSE)]
  
  ggs = read.xlsx(paste0('../data/gene_list_for_TM3Seq.xlsx'), sheet = 1, colNames = FALSE)
  colnames(ggs) = c('pathway', 'gene')
  ggs = ggs[, c(2, 1)]
  
  p = ggs$pathway[1]
  for(n in 2:nrow(ggs))
  {
    if(is.na(ggs$pathway[n])) {
      ggs$pathway[n] = p  
    }else{
      if(ggs$pathway[n] != p) p = ggs$pathway[n]
    }
  }
  
  Add.GoTerm.gene.list = TRUE
  if(Add.GoTerm.gene.list){
    for(n in 1:length(xx.files))
    {
      # n = 1
      xx = read.xlsx(xx.files[n])
      xx = unique(xx$Symbol)
      xx = cbind(xx, rep(unlist(strsplit(as.character(gsub('.xlsx', '', basename(xx.files[n]))), '_'))[4], length(xx)))
      colnames(xx) = c('gene', 'pathway')
      ggs = rbind(ggs, xx)
    }
  }
  
  ggs = rbind(c('Id3', 'BMP'), ggs)
  ggs = rbind(c('Id1', 'BMP'), ggs)
  ggs = rbind(c('Smad6', 'BMP'), ggs)
  ggs = rbind(c('Foxa2', NA), ggs)
  
  gg.unique = unique(ggs$gene)
  ggs = ggs[match(gg.unique, ggs$gene), ]
  
  xx = read.xlsx('../data/pbio.2002117.s012.xlsx')
  xx = xx[, c(2,3)]
  xx = xx[!is.na(xx$geneName), ]
  colnames(xx) = c('gene', 'pathway')
  
  firstup <- function(x) {
    x <- tolower(x)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }
  
  xx$gene = sapply(xx$gene, firstup)
  ggs = rbind(ggs, xx)
  gg.unique = unique(ggs$gene)
  ggs = ggs[match(gg.unique, ggs$gene), ]
  ggs = ggs[order(ggs$pathway), ]
  
  saveRDS(ggs, file = paste0(RdataDir, '/TM3_examplesGenes_withGOterm.rds'))
  
  #write.csv(ggs, file = paste0(resDir, '/genes_signalingPathways.csv'), row.names = FALSE)
}


########################################################
########################################################
# Section : test pooling of positive and negative cells
# 
########################################################
########################################################
Test.pooling.TM3.positive.negative.cells = function()
{
  Test.pooling.negative.positive.cells = FALSE
  if(Test.pooling.negative.positive.cells){
    fpm = fpm(dds, robust = TRUE)
    
    rrs = readRDS(file = paste0(RdataDir, 'facs_positive_negative_ratios.rds'))
    rrs$Tube.Name = as.character(rrs$Tube.Name)
    rrs$Tube.Name = gsub('[-]', '_', rrs$Tube.Name)
    rrs$Tube.Name = gsub('Fgf', 'FGF', rrs$Tube.Name)
    rrs$Tube.Name = gsub('PD03', 'PD', rrs$Tube.Name)
    rrs$Tube.Name = gsub('Chiron', 'CHIR', rrs$Tube.Name)
    
    pools = matrix(NA, nrow = nrow(fpm), ncol = length(unique(design.matrix$condition.rep)))
    colnames(pools) = unique(design.matrix$condition.replicate)
    rownames(pools) = rownames(fpm)
    
    for(n in 1:ncol(pools))
    {
      # n = 1
      ratio = rrs$GFP_pos..Parent[which(rrs$Tube.Name == colnames(pools)[n])]/100  
      
      jj1 = which(design.matrix$condition.replicate == colnames(pools)[n] & design.matrix$cells == 'Foxa2.pos')
      jj2 = which(design.matrix$condition.replicate == colnames(pools)[n] & design.matrix$cells == 'Foxa2.neg')
      
      if(length(jj1) == 1 & length(jj2) == 1){
        cat(n, ' - ', colnames(pools)[n], '- positive cell ratio ', ratio,  ' with column ', jj1, jj2, '\n')
        pools[,n] = ratio * fpm[,jj1] + (1 - ratio) * fpm[, jj2]
      }
      
      
    }
    
    pools = pools[, which(colnames(pools) != 'LDN_1')]
    
    cc.pools = colnames(pools)
    cc.pools = sapply(cc.pools, function(x) unlist(strsplit(as.character(x), '_'))[1])
    
    library(factoextra)
    ntop = 500
    xx = as.matrix(log2(pools + 2^-4))
    vars = apply(xx, 1, var)
    xx = xx[order(-vars), ]
    xx = xx[1:ntop, ]
    
    res.pca <- prcomp(t(xx), scale = TRUE)
    #res.var <- get_pca_var(res.pca)
    
    fviz_pca_ind(res.pca,
                 col.ind = "cos2", # Color by the quality of representation
                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                 repel = TRUE     # Avoid text overlapping
    )
    
    save(pools, cc.pools, file = paste0(RdataDir, '/TM3_pooled.positive.negative_', Counts.to.Use, version.analysis, '.Rdata'))
    
    g = 'Foxa2'
    j = grep(g, rownames(pools))
    
    log2(fpm[j, grep('neg', colnames(fpm))] + 2^-6)
    log2(fpm[j, grep('pos', colnames(fpm))] + 2^-6)
    log2(pools[j, ] + 2^-6)
    
    
  }
  
  ##########################################
  # check the ratios between posive and negative cells  
  ##########################################
  design.matrix$condition.rep = paste0(design.matrix$condition, '_', design.matrix$Triplicate.No.)
  ratios = matrix(NA, nrow = nrow(fpm), ncol = length(unique(design.matrix$condition.rep)))
  colnames(ratios) = unique(design.matrix$condition.rep)
  rownames(ratios) = rownames(fpm)
  
  for(n in 1:nrow(ratios))
  {
    # n = which(rownames(ratios) == 'Chrd')
    cat(n, '\n')
    for(j in 1:ncol(ratios))
    {
      k.pos = which(design.matrix$condition.rep == colnames(ratios)[j] & design.matrix$cells == 'Foxa2.pos')
      k.neg = which(design.matrix$condition.rep == colnames(ratios)[j] & design.matrix$cells == 'Foxa2.neg')
      ratios[n, j] = fpm[n, k.pos] - fpm[n, k.neg]
    }
    
  }
  
}