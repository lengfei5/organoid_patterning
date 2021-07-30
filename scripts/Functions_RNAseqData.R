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
  
  ##########################################
  # manually add some examples 
  ##########################################
  spg = readRDS(file = '../results/Rdata/TM3_examplesGenes_withGOterm.rds')
  spg = rbind(c('Olig2', NA), spg)
  spg = rbind(spg, c('Shh', 'SHH'))
  spg = rbind(spg, c('Runx1', 'BMP'))
  
  spg = spg[match(unique(spg$gene), spg$gene), ]
  saveRDS(spg, file = paste0('../results/Rdata/curated_signaling.pathways_gene.list.rds'))
  
  ##########################################
  # add SHH pathways  
  ##########################################
  genes = readRDS(file = '../results/Rdata/curated_signaling.pathways_gene.list.rds')
  
  xx = read.xlsx('../data/GO_term_summary_SHH.xlsx')
  xx = unique(xx$Symbol)
  xx = cbind(xx, rep('SHH', length(xx)))
  colnames(xx) = colnames(genes)
  genes = rbind(genes, xx)  
  #genes = genes[, genes$gene), ]
  
  saveRDS(genes, file = paste0('../results/Rdata/curated_signaling.pathways_gene.list_v2.rds'))
  #write.csv(ggs, file = paste0(resDir, '/genes_signalingPathways.csv'), row.names = FALSE)
}

Save.median.transcript.lengths.for.each.gene = function()
{
  annot = read.delim(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/', 
                            'ens_BioMart_GRCm38.p6.txt'), sep = '\t', header = TRUE)
  
  annot = annot[which(annot$Gene.type == 'protein_coding'), ]
  annot = annot[!is.na(match(annot$Chromosome.scaffold.name, as.character(c(1:19, 'X', 'Y')))), ]
  
  gg.uniq = unique(annot$Gene.stable.ID)
  mm = match(gg.uniq, annot$Gene.stable.ID)
  annots = annot
  
  genes = data.frame(gg.uniq, annots$Gene.name[mm], gene.types = annots$Gene.type[mm], stringsAsFactors = FALSE)
  genes = genes[which(!is.na(genes[, 3])), ]
  
  #genes = genes[match(unique(genes[,3]), genes[,3]), ]
  colnames(genes) = c('ensID', 'gene', 'genetype')
  genes$length = annot$Transcript.length..including.UTRs.and.CDS.[match(genes$ensID, annots$Gene.stable.ID)]
  genes$nb.transcript = annot$Transcript.count[match(genes$ensID, annot$Gene.stable.ID)]
  
  for(n in 1:nrow(genes))
    #for(n in 1:200)
  {
    # n = 1
    if(n%%100 ==0) cat(n, '\n')
    if(genes$nb.transcript[n] > 1){
      genes$length[n] = median(annots$Transcript.length..including.UTRs.and.CDS.[which(annots$Gene.stable.ID == genes$ensID[n])])
    }
  }
  
  saveRDS(genes, file = paste0(RdataDir, 'mm10_ens_BioMart_GRCm38.p6_ensID.geneSymbol.length.rds'))
  
}

########################################################
########################################################
# Section : test pooling of positive and negative cells
# two methods for pooling, with counts or normalized fpm
########################################################
########################################################
Test.pooling.TM3.positive.negative.cells = function(dds, pooling = 'counts', design.matrix = design.matrix)
{
  # Import the positive and negative ratios
  rrs = readRDS(file = paste0(RdataDir, 'facs_positive_negative_ratios.rds'))
  rrs$Tube.Name = as.character(rrs$Tube.Name)
  rrs$Tube.Name = gsub('[-]', '_', rrs$Tube.Name)
  rrs$Tube.Name = gsub('Fgf', 'FGF', rrs$Tube.Name)
  rrs$Tube.Name = gsub('PD03', 'PD', rrs$Tube.Name)
  rrs$Tube.Name = gsub('Chiron', 'CHIR', rrs$Tube.Name)
  
  pools = matrix(NA, nrow = nrow(dds), ncol = length(unique(design.matrix$condition.rep)))
  colnames(pools) = unique(design.matrix$condition.replicate)
  rownames(pools) = rownames(dds)
  
  if(pooling == 'counts' ) {sorted = counts(dds)}
  if(pooling == 'fpm'){ sorted = fpm(dds)}
  
  for(n in 1:ncol(pools))
  {
    # n = 1
    ratio = rrs$GFP_pos..Parent[which(rrs$Tube.Name == colnames(pools)[n])]/100  
    
    jj1 = which(design.matrix$condition.replicate == colnames(pools)[n] & design.matrix$cells == 'Foxa2.pos')
    jj2 = which(design.matrix$condition.replicate == colnames(pools)[n] & design.matrix$cells == 'Foxa2.neg')
    
    if(length(jj1) == 1 & length(jj2) == 1){
      cat(n, ' - ', colnames(pools)[n], '- positive cell ratio ', ratio,  ' with column ', jj1, jj2, '\n')
      pools[,n] = ratio * sorted[,jj1] + (1 - ratio) * sorted[, jj2]
    }
    
  }
  
  pools = pools[, which(colnames(pools) != 'LDN_1')]
  
  cc.pools = colnames(pools)
  cc.pools = sapply(cc.pools, function(x) unlist(strsplit(as.character(x), '_'))[1])
  
  if(pooling == 'counts'){
    pools = ceiling(as.matrix(pools))
    ddx <- DESeqDataSetFromMatrix(pools, DataFrame(condition = cc.pools), design = ~ condition)
    
    ss = rowSums(counts(ddx))
    
    hist(log2(ss), breaks = 200, main = 'log2(sum of reads for each gene)')
    # 
    length(which(ss > quantile(ss, probs = 0.5)))
    dd0 = ddx[ss > quantile(ss, probs = 0.6) , ]
    dd0 = estimateSizeFactors(dd0)
    sizefactors.UQ = sizeFactors(dd0)
    # 
    # plot(sizeFactors(dd0), colSums(counts(dds)), log = 'xy')
    # text(sizeFactors(dd0), colSums(counts(dds)), colnames(dd0), cex =0.4)
    # 
    # use UQ normalization from edgeR
    #library(edgeR)
    #dge2 <- DGEList(raw)
    #dge2 <- calcNormFactors(dge2, method = "upperquartile")
    #dge2$samples
    #sizefactors.UQ = as.data.frame(dge2$samples) 
    #sizefactors.UQ = sizefactors.UQ$lib.size * sizefactors.UQ$norm.factors/median(sizefactors.UQ$lib.size)
    
    # cutoff.gene = 50
    # cat(length(which(ss > cutoff.gene)), 'genes selected \n')
    # 
    # dds <- dds[ss > cutoff.gene, ]
    
    # normalization and dimensionality reduction
    sizeFactors(ddx) = sizefactors.UQ
    cpm = fpm(ddx, robust = TRUE)
    
    #save(fpm, design, file = paste0(tfDir, '/RNAseq_fpm_fitered.cutoff.', cutoff.gene, '.Rdata'))
    vsd <- varianceStabilizingTransformation(ddx, blind = FALSE)
    
    pca=plotPCA(vsd, intgroup = c('condition'), returnData = TRUE, ntop = 2000)
    pca2save = as.data.frame(pca)
    
    ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition))  + 
      geom_point(size=2.5) + 
      #geom_text(hjust = 0.4, nudge_y = 0.7, size=3) + 
      geom_text_repel()
      #geom_label_repel()
    plot(ggp)
    
    library(factoextra)
    ntop = 500
    xx = as.matrix(log2(cpm + 2^-4))
    means = apply(xx, 1, mean)
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
    
    save(ddx, cc.pools, file = paste0(RdataDir, '/TM3_pooled.pos.neg_ddx_cc.pools_', Counts.to.Use, version.analysis, '.Rdata'))
    
  }
  
  if(pooling == 'fpm'){
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

##########################################
# Compare the TM3 data and Elena's RNA-seq data 
##########################################
Compare.TM3.and.RNAseq.timeSeries.sortedDay5 = function()
{
  fpm = fpm(dds, robust = TRUE)
  cpm = log2(fpm + 2^-4)
  cc = paste0(design.matrix$condition, '_', design.matrix$cells)
  
  # load samples with pooled negative and positive cells
  #load(file = paste0(RdataDir, '/TM3_positive.negative.pooled_', Counts.to.Use, version.analysis, '.Rdata'))
  load(file = paste0(RdataDir, '/TM3_pooled.positive.negative_', Counts.to.Use, version.analysis, '.Rdata'))
  cc.pools = paste0(cc.pools, '_pooled')
  pools = log2(pools + 2^-4)
  colnames(pools) = paste0(colnames(pools), '.pooled')
  
  cpm = cbind(cpm, pools)
  cc = c(cc, cc.pools)
  
  #load(file = paste0('../results/Rdata/RNAseq_timeSeries_sortedDay5_count_design_geneSymbol_dds.Rdata'))
  #design0 = design
  #dds0 = dds
  #save(dds0, design0, file = paste0('../results/Rdata/RNAseq_timeSeries_sortedDay5_count_design_geneSymbol_dds_backup.Rdata'))
  load(file = paste0('../results/Rdata/RNAseq_timeSeries_sortedDay5_count_design_geneSymbol_dds_backup.Rdata')) 
  fpm0 = fpm(dds0)
  
  kk = c(grep('RA', design0$condition), which(design0$condition=='12h'))
  #kk = c(1:nrow(design0))
  
  cc0 = as.character(design0$condition[kk])
  xx = log2(fpm0[, kk] + 2^-4) 
  
  bc = c(rep(1, length(cc)), rep(0, length(cc0)))
  cc[which(cc == 'RA_pooled')] = '12h_RA'
  cc0[which(cc0 == '12h')] = 'N2B27_pooled'
  cc = c(cc, as.character(cc0))
  mm = match(rownames(cpm), rownames(xx))
  
  cpm = cbind(cpm[which(!is.na(mm)), ], xx[mm[which(!is.na(mm))], ])
  
  jj = c(grep('N2B27|pos|neg|pooled', cc, invert = TRUE), grep("RA_Foxa2.pos|RA_Foxa2.neg|N2B27_pooled", cc))
  
  cc = cc[jj]
  cpm = cpm[,jj]
  bc = bc[jj]
  
  require("sva")
  bc = as.factor(bc)
  mod = model.matrix(~ as.factor(conds), data = data.frame(conds = cc))
  
  yy = ComBat(dat=cpm, batch=bc, mod=mod, par.prior=TRUE, ref.batch = 0)    
  
  #yy = cpm
  ntop = 5000
  xx = as.matrix(yy)
  vars = apply(xx, 1, var)
  xx = xx[order(-vars), ]
  xx = xx[1:ntop, ]
  #par(mfrow = c(2,2))
  #pairs(xx[, c(1:4)])
  #plot.pair.comparison.plot(xx[, c(1:4)], linear.scale = FALSE)
  #plot.pair.comparison.plot(xx[, c(9:16)], linear.scale = FALSE)
  
  res.pca <- prcomp(t(xx), scale = TRUE)
  pcas = data.frame(res.pca$x[, c(1:2)], condition = cc)
  pcas = data.frame(pcas, name = rownames(pcas))
  
  ggplot(data=pcas, 
         aes(PC1, PC2, color= condition, label = name))  + 
    geom_point(size=3) + 
    geom_text(hjust = 1, nudge_y = 0.5, size=3) + 
    ggsave(paste0(resDir, "/PCAplot_ElenaRNAseq_pooledRA.pos.neg.pdf"), width=18, height = 12)
  
  #plot(res.pca$x[, c(1,2)])
  # fviz_pca_ind(res.pca,
  #              col.ind = "cos2", # Color by the quality of representation
  #              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  #              repel = TRUE     # Avoid text overlapping
  #)
  
  # comapre the foxa2 positive and negative cells for RA in day3 and day5
  xx1 = apply(cpm[, grep('s48h_RA_GFPp', colnames(cpm))], 1, median)- apply(cpm[ ,grep('s48h_RA_AF', colnames(cpm))], 1, median) 
  xx2 = apply(cpm[, grep('RA_[[:digit:]]_Foxa2.pos', colnames(cpm))], 1, median) - 
    apply(cpm[, grep('RA_[[:digit:]]_Foxa2.neg', colnames(cpm))], 1, median)
  
}


Compare.TM3.sortedDay3.RNAseq.sortedDay5 = function(dds)
{
  kk = which(dds$conds == 'RA_Foxa2.pos' | dds$conds == 'RA_Foxa2.neg')
  ddx = dds[,kk]
  ddx$conds <- droplevels(ddx$conds)
  
  load(file = paste0('../results/Rdata/RNAseq_Foxa.positive_vs_neg.Day5.Rdata'))
  
  gene.sels = unique(intersect(rownames(ddx), rownames(dds1)))
  ddx = ddx[match(gene.sels, rownames(ddx)), ]
  dds1 = dds1[match(gene.sels, rownames(dds1))]
  
  raw = cbind(counts(ddx), counts(dds1))
  cc = c(as.character(ddx$conds), as.character(dds1$condition))
  
  ddx <- DESeqDataSetFromMatrix(raw, DataFrame(condition = cc), design = ~ condition)
  ddx = estimateSizeFactors(ddx)
  vsd <- varianceStabilizingTransformation(ddx, blind = FALSE)
  
  pca=plotPCA(vsd, intgroup = c('condition'), returnData = TRUE, ntop = 500)
  pca2save = as.data.frame(pca)
  ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition))  + 
    geom_point(size=4) + 
    geom_text(hjust = 0.2, nudge_y = 0.25, size=3)
  
  plot(ggp)
  
}

########################################################
# saturation curve from rseqc
# the r code from rseqc output
########################################################
sequencing.saturation.analysis = function()
{
  Sequence.saturation.analysis = FALSE
  if(Sequence.saturation.analysis){
    rseqc.file = list.files('../Data/R10724_rnaseq/saturation_rseqc', pattern = 'junctionSaturation_plot.r', 
                            full.names = TRUE)
    library(stringr)
    
    yy = c()
    for(n in 1:length(rseqc.file))
    {
      cat(n, '\n')
      xx = read.delim(rseqc.file[n])
      xx = xx[grep('y=', xx[, 1 ]), ]
      #xx = gsub('y=c', '', xx)
      # Get the parenthesis and what is inside
      k <- str_extract_all(xx, "\\([^()]+\\)")[[1]]
      # Remove parenthesis
      k <- substring(k, 2, nchar(k)-1)
      #k = gsub('["]', '', k)
      k = as.numeric(unlist(strsplit(as.character(k), ',')))
      yy = rbind(yy, k)
    }
    
    rownames(yy) = gsub('_junction.junctionSaturation_plot.r', '', basename(rseqc.file))
    
    
    pdfname = paste0(resDir, '/saturation_curve_rseqc_knownJunctions.pdf')
    pdf(pdfname, width = 16, height = 8)
    par(cex = 1.0, las = 1,  mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0), tcl = -0.3)
    
    yy = yy[which(rownames(yy) != '136150_TTAACCTTCGAGGCCAGACA_HNF3KDSXY_3_20201223B_20201223'), ]
    span = 0.75
    # saturation curve with nb of peaks
    xlims = c(0, 120)
    ylims = range(yy/10^3)
    frac = c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100)/100
    
    library(RColorBrewer)
    cols = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(nrow(yy))
    plot(0, 0, xlim = xlims, ylim = ylims, type ='n', xlab = 'nb of TOTAL reads (Million)', 
         ylab = 'nb of known junctions (K)', main = paste0('saturation curve from rseqc'))
    abline(v = c(20, 30, 40, 50), col = 'blue', lwd = 1.0, lty =2)
    
    #legend('topleft', legend = sample.uniq, col = cols, bty = 'n', lwd = 2.0, cex = 0.7)
    
    for(n in 1:nrow(yy))
    {
      # n = 1
      cat(n, '\n')
      
      kk = which(design$sample == rownames(yy)[n])
      
      satt = data.frame(nb.reads = design$total.reads[kk]*frac/10^6, nb.junctions = yy[n, ]/10^3)
      
      points(satt[,1], satt[,2], type= 'p', col = cols[n])
      loessMod <- loess(nb.junctions ~ nb.reads, data=satt, span=span)
      smoothed <- predict(loessMod)
      lines(smoothed, x=satt$nb.reads, col=cols[n], lwd = 3.0)
      
      text(satt[nrow(satt), 1], smoothed[length(smoothed)], labels = paste0(design$fileName[kk], '_', design$sampleID[kk]), 
           cex = 0.7, pos = 4, offset = 0.2)
      
    }
    
    dev.off()
    
  }
  
}
