##########################################################################
##########################################################################
# Project:
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Sat Jun 19 16:34:59 2021
##########################################################################
##########################################################################
########################################################
########################################################
# Section I : data processing and sequencing quality controls
# 
########################################################
########################################################
rm(list = ls())
RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

# setup for data import and sequencing QCs
version.analysis = '_R11601_TM3_20210619'

resDir = paste0("../results/RNAseq_perturbation", version.analysis)
RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/patterning_organoid/R11601/'

Collect.QCs.stat = TRUE
Counts.to.Use = "UMI"

##########################################
# prepare design matrix, count table and QC table
##########################################
design = read.csv(paste0(dataDir, 'sampleInfo.csv'), header = TRUE)
design$condition = as.character(design$condition)
design$cells = design$condition

design$cells[grep('[+]', design$cells)] = 'Foxa2+'
design$cells[grep('[-]', design$cells)] = 'Foxa2-'
design$condition = sapply(design$condition, function(x) unlist(strsplit(as.character(x), '[_]'))[1])

if(Collect.QCs.stat){
  MultiQC.Dir = paste0(dataDir, 'multiqc_data/')
  
  stats = read.delim(paste0(MultiQC.Dir, 'multiqc_general_stats.txt'), sep = '\t', header = TRUE)
  alignment = read.delim(paste0(MultiQC.Dir, 'multiqc_star.txt'), sep = '\t')
  
  ##########################################
  # process the stat table
  ##########################################
  keep = data.frame(matrix(NA, nrow = nrow(design), ncol = 11))
  colnames(keep) = c('sampleID', 'pct.duplication', 'pct.GC', 'avg.seq.length', 'total.reads', 
                     'pct.assign', 'assigned.reads', 'alignment.rate', 'trimmed.reads', 'unique.aligned', 'multimapper')
  keep$sampleID = design$sampleID
  
  stats = data.frame(stats, stringsAsFactors = FALSE)
  alignment = data.frame(alignment, stringsAsFactors = FALSE)
  
  
  for(n in 1:nrow(keep))
  {
    jj = grep(keep$sampleID[n], stats$Sample)
    jj = jj[grep('R1', stats$Sample[jj])]
    jj1 = jj[grep('R1_trimmed_sorted_umiDedup', stats$Sample[jj], invert = TRUE)]
    jj2 = jj[grep('R1_trimmed_sorted_umiDedup', stats$Sample[jj])]
    keep[n, c(2, 3, 4, 5)] = stats[jj1, c(2, 3, 4, 6)]
    keep[n, c(6, 7)] = stats[jj2, c(9, 10)]
    jj3 = grep(keep$sampleID[n], alignment$Sample)
    keep$alignment.rate[n] = alignment$uniquely_mapped_percent[jj3]
    keep$trimmed.reads[n] = alignment$total_reads[jj3]
    keep$unique.aligned[n] = alignment$uniquely_mapped[jj3]
    keep$multimapper[n] = alignment$multimapped[jj3]
     
  }
  
  #stats = stats[, c(1, 2, 3, 4, 6, 8, 9, 10)]
  #stats = data.frame(stats, alignment[match(stats$sample, alignment$Sample), c(2, 4, 5)], stringsAsFactors = FALSE)
  #colnames(stats)[c(9:11)] = c()
  #stats = stats[, c(1:5, 9, 8, 10, 11, 7)]
  
  #design = design[order(design$fileName), ]
  # ii = c()
  # jj = c()
  # for(n in 1:nrow(design))
  # {
  #   # n = 1;
  #   cat(n, '\n')
  #   kk = grep(design$sampleID[n], stats$sample)
  #   ii = c(ii, rep(n, length(kk)))
  #   jj = c(jj, kk)
  #   #kk = c(kk, grep(design$sampleID[n], stats$sample))
  #   #kk = c(kk, which(design$sampleID == ))
  # }
  
  xx = data.frame(design, keep[, -1], stringsAsFactors = FALSE)
  #xx = xx[order(xx$fileName), ]
  
  write.csv(xx, file = paste0(resDir, '/sampleInfos_QCs.stats.csv'), row.names = FALSE)
  
  design = xx
  
  saveRDS(design, file = paste0(RdataDir, 'sampleInfo_QC.stats.rds'))
    
}

########################################################
# saturation curve from rseqc
# the r code from rseqc output
########################################################
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

##################################################
## Import UMI count table
##################################################
source(RNA.functions)
source(RNA.QC.functions)

design = readRDS(file = paste0(RdataDir, 'sampleInfo_QC.stats.rds'))

colnames(design)[1] = 'SampleID'
design$conds = paste0(design$condition, '_', design$cells)

# prepare the data table 
Dir_umi = paste0(dataDir, "htseq_counts_BAMs_umi")
Dir_read = paste0(dataDir, "htseq_counts_BAMs")

aa1 <- list.files(path = Dir_umi, pattern = "*umiDedup.txt", full.names = TRUE)
aa1 = merge.countTables.htseq(aa1)
colnames(aa1)[-1] = paste0(colnames(aa1)[-1], ".UMI")

aa2 <- list.files(path = Dir_read, pattern = "*.txt", full.names = TRUE)
aa2 = merge.countTables.htseq(aa2)
colnames(aa2)[-1] = paste0(colnames(aa2)[-1], ".readCount")

aa <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), list(aa1, aa2))

## compare read counts vs. umi counts
source(RNA.functions)
Compare.UMI.vs.readCounts = TRUE
if(Compare.UMI.vs.readCounts){
    pdfname = paste0(resDir, "/readCounts_vs_UMI_normalized", version.analysis, ".pdf")
    pdf(pdfname, width = 10, height = 8)

    compare.readCount.UMI(design[, c(1:3)], aa, normalized = FALSE)

    dev.off()
}


if(Counts.to.Use == 'readCounts'){
   all = process.countTable(all=aa, design = design, special.column = ".readCount", ensToGeneSymbol = TRUE)
 }else{
   if(Counts.to.Use == "UMI"){
     all = process.countTable(all=aa, design = design, special.column = "UMI", ensToGeneSymbol = TRUE)
   }else{
     cat("Error : no counts found for ", Counts.to.Use, "for miRNAs \n")
   }
 }


all = process.countTable(all=all, design = design[,c(1, 14)], special.column = NULL, ensToGeneSymbol = FALSE)

all = all[grep('^__', all$gene, invert = TRUE), ]

save(design, all, file=paste0(RdataDir, '/Design_Raw_readCounts', Counts.to.Use,  version.analysis, '.Rdata'))

########################################################
########################################################
# Section II : analyze the RNA-seq data 
# 
########################################################
########################################################

##########################################
# gene names converted from ensID to gene symbol 
##########################################
load(file=paste0(RdataDir, '/Design_Raw_readCounts', Counts.to.Use,  version.analysis, '.Rdata'))

design$cells[grep('[+]', design$cells)] = 'Foxa2.pos'
design$cells[grep('[-]', design$cells)] = 'Foxa2.neg'
design$conds = paste0(design$condition, '_', design$cells)

colnames(all)[-1] = paste0(design$condition, '_', design$cells, '_', design$SampleID)

annot = read.delim(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/', 
                       'ens_BioMart_GRCm38.p6.txt'), sep = '\t', header = TRUE)

Select.proteinCoding.genes = TRUE
if(Select.proteinCoding.genes){
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
    
    all = all[!is.na(match(all$gene, annot$Gene.stable.ID)), ]
    
    rownames(all) = all$gene
    all$gene = annot$Gene.name[match(all$gene, annot$Gene.stable.ID)]
    
    all = all[!is.na(all$gene), ]
    gg.uniq = unique(all$gene)
    all = all[match(gg.uniq, all$gene), ]
    rownames(all) = all$gene
    all = all[, -1]
}


Add.more.sample.details = TRUE
if(Add.more.sample.details){
  library("openxlsx")
  xx = read.xlsx(paste0(dataDir, 'Detailed_sample_information_with_SeqIDs.xlsx'), sheet = 1)
  xx = data.frame(design, xx[match(design$SampleID, xx$Seq.ID), ], stringsAsFactors = FALSE)
  design = xx
  
  colnames(design)[grep('Triplicate', colnames(design))] = 'triplicate.nb'
  
  xx = read.csv(paste0(dataDir, '20210429_Venus_F02-Batch_Analysis_29042021122425.csv'))
  saveRDS(xx, file = paste0(RdataDir, 'facs_positive_negative_ratios.rds'))
  
}

save(design, all, file = paste0(RdataDir, '/design.detailed_RawUMI_', Counts.to.Use, version.analysis, '.Rdata'))

##########################################
# QCs of replicates and conditions
##########################################
load(file = paste0(RdataDir, '/design.detailed_RawUMI_', Counts.to.Use, version.analysis, '.Rdata'))

QC.for.cpm = FALSE
if(QC.for.cpm){
  
  source(RNA.functions)
  source(RNA.QC.functions)
    
  raw = as.matrix(all)
  
  kk = which(design$SampleID != '161040' & design$condition != 'N2B27')
  #raw = raw[, -kk]
  
  ss = apply(as.matrix(raw), 1, sum)
  raw = raw[which(ss >100), ]
  
  pdfname = paste0(resDir, "/Data_qulity_assessment_umiThreshold.100_noNegativeControls", version.analysis, ".pdf")
  pdf(pdfname, width = 12, height = 10)
  
  Check.RNAseq.Quality(read.count=raw[, kk], design.matrix = design[kk, c(1, 2, 3)], 
                       lowlyExpressed.readCount.threshold=500)
  
  dev.off()
  
}


########################################################
########################################################
# Section : Analysis with DESeq2
# 
########################################################
########################################################
require(ggplot2)
require(DESeq2)
require(gridExtra)
library(dplyr)
library(patchwork)

Counts.to.Use = 'UMI'

load(file = paste0(RdataDir, '/design.detailed_RawUMI_', Counts.to.Use, version.analysis, '.Rdata'))
genes = readRDS(file = paste0(RdataDir, 'mm10_ens_BioMart_GRCm38.p6_ensID.geneSymbol.length.rds'))
ggs = readRDS(file = paste0(RdataDir, '/TM3_examplesGenes_withGOterm.rds'))

sels = which(design$SampleID != '161040')

design.matrix = design[sels, ]
design = design[sels, ]
rm(design)

raw = all[, sels]
design.matrix$condition.replicate = paste0(design.matrix$condition, '_', design.matrix$triplicate.nb)
colnames(raw) = paste0(design.matrix$condition.replicate, '_', design.matrix$cells, '_', design.matrix$SampleID)

dds <- DESeqDataSetFromMatrix(raw, DataFrame(design.matrix), design = ~ conds)


ss = rowSums(counts(dds))

hist(log2(ss), breaks = 200, main = 'log2(sum of reads for each gene)')

dd0 = dds[ss > quantile(ss, probs = 0.6) , ]
dd0 = estimateSizeFactors(dd0)
sizefactors.UQ = sizeFactors(dd0)

plot(sizeFactors(dd0), colSums(counts(dds)), log = 'xy')
text(sizeFactors(dd0), colSums(counts(dds)), colnames(dd0), cex =0.4)

# use UQ normalization from edgeR
#library(edgeR)
#dge2 <- DGEList(raw)
#dge2 <- calcNormFactors(dge2, method = "upperquartile")
#dge2$samples
#sizefactors.UQ = as.data.frame(dge2$samples) 
#sizefactors.UQ = sizefactors.UQ$lib.size * sizefactors.UQ$norm.factors/median(sizefactors.UQ$lib.size)

cutoff.gene = 50
cat(length(which(ss > cutoff.gene)), 'genes selected \n')

dds <- dds[ss > cutoff.gene, ]

# normalization and dimensionality reduction
sizeFactors(dds) = sizefactors.UQ
fpm = fpm(dds, robust = TRUE)

#save(fpm, design, file = paste0(tfDir, '/RNAseq_fpm_fitered.cutoff.', cutoff.gene, '.Rdata'))
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('condition', 'cells'), returnData = TRUE, ntop = 500)
#print(pca)d
pca2save = as.data.frame(pca)
ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition, shape = cells))  + 
  geom_point(size=4) + 
  geom_text(hjust = 0.2, nudge_y = 0.25, size=3)

plot(ggp) + ggsave(paste0(resDir, "/PCAplot_withControls_UQ.norm_ntop500.pdf"), width=18, height = 12)

ggp = ggplot(data=pca2save[grep('N2B27', pca2save$condition, invert = TRUE), ], 
             aes(PC1, PC2, label = name, color= condition, shape = cells))  + 
  geom_point(size=4) + 
  geom_text(hjust = 0.2, nudge_y = 0.25, size=3)

plot(ggp) + ggsave(paste0(resDir, "/PCAplot_withoutControls_UQ.norm_ntop500.pdf"), width=18, height = 12)

##########################################
# try to put the TM3 data in the context of time series data
# test the relationship between perturbation of positive and negative cells, pooled cells, and time seires and day 5 sorted cells 
# 
##########################################
Compare.with.pooledCells.timeSeries.sortedDay5 = FALSE
if(Compare.with.pooledCells.timeSeries.sortedDay5){
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

########################################################
########################################################
# Section : Pairwise comparisons 
# Here we first identify significantly different genes in RA positive and negative genes
# RA positive vs negative
# perturbed-positive vs RA-positive; perturbed-negative vs RA-negative
########################################################
########################################################
Calculate.pairwise.comparisons = FALSE
if(Calculate.pairwise.comparisons){
  
  Make.pairwise.comparisons = FALSE
  if(Make.pairwise.comparisons){
    library(org.Mm.eg.db)
    library(enrichplot)
    library(clusterProfiler)
    library(ggplot2)
    library(stringr)
    
    cc = unique(design.matrix$condition)
    cc = cc[which(cc != "N2B27")]
    
    bgs = rownames(dds)
    bgs.df <- bitr(bgs, fromType = "SYMBOL",
                   toType = c("ENSEMBL", "ENTREZID"),
                   OrgDb = org.Mm.eg.db)
    
    pdfname = paste0(resDir, '/TM3_pairwiseComparisons_enrichGO_pathways.pdf')
    pdf(pdfname,  width = 18, height = 16)
    #par(cex = 1.0, las = 1, mgp = c(3,2,0), mar = c(6,6,2,0.2), tcl = -0.3)
    
    for(n in 1:length(cc))
    #for(n in 1:2)
    {
      # n = 3
      cat(as.character(cc[n]), '\n')
      if(cc[n] == 'RA'){
        kk = which(design.matrix$condition == cc[n])
      }else{
        kk = which(design.matrix$condition == cc[n] | design.matrix$condition == 'RA')
      }
      dds1 = dds[, kk]
      dds1$conds <- droplevels(dds1$conds)
      dds1 <- estimateDispersions(dds1, fitType = 'parametric')
      plotDispEsts(dds1, ymin = 10^-3, main = cc[n])
      #dds = estimateDispersions(dds)
      #plotDispEsts(dds)
      dds1 <- nbinomWaldTest(dds1)
      resultsNames(dds1)  
      
      if(cc[n] == 'RA'){
        res1 = results(dds1, contrast=c("conds", 'RA_Foxa2.pos', 'RA_Foxa2.neg'), alpha = 0.05)
        res1 <- lfcShrink(dds1, coef=2, type="normal")
        res1 = as.data.frame(res1)
        gg.signif = rownames(res1)[which(res1$padj < 0.05)]
        
        colnames(res1) = paste0(colnames(res1), '_RA.pos_vs_RA.neg')
        
        # check the siganificant ones and overlapping with SP genes
        #xx = data.frame(res[order(res1$pvalue), ])
        #xx = xx[which(xx$padj < 0.05), ]
        #xx = xx[!is.na(match(rownames(xx), ggs$gene)), ]
        
      }else{
        res1 = results(dds1, contrast=c("conds", paste0(cc[n], '_Foxa2.pos'), 'RA_Foxa2.pos'), alpha = 0.05)
        res1 <- lfcShrink(dds1, contrast=c("conds", paste0(cc[n], '_Foxa2.pos'), 'RA_Foxa2.pos'), type="normal")
        res1 = as.data.frame(res1)
        
        res2 = results(dds1, contrast=c("conds", paste0(cc[n], '_Foxa2.neg'), 'RA_Foxa2.neg'), alpha = 0.05)
        res2 <- lfcShrink(dds1, type="normal", contrast=c("conds", paste0(cc[n], '_Foxa2.neg'), 'RA_Foxa2.neg'))
        res2 = as.data.frame(res2)
        
        gg.signif = rownames(res1)[which(res1$padj < 0.1 | res2$padj < 0.1)]
        cat(length(gg.signif), ' significant genes for enrichGO \n')
        
        gene.df <- bitr(gg.signif, fromType = "SYMBOL", toType = c("ENSEMBL", "ENTREZID"), OrgDb = org.Mm.eg.db)
        ego <-  enrichGO(gene         = gene.df$ENSEMBL,
                         universe     = bgs.df$ENSEMBL,
                         #OrgDb         = org.Hs.eg.db,
                         OrgDb         = org.Mm.eg.db,
                         keyType       = 'ENSEMBL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05, readable = TRUE)
        
        p0 = barplot(ego, showCategory=50, main = cc[n]) 
        plot(p0)
        
        ii = grep('pathway', ego[, 2])
        #xx = ego[ii, ]
        write.table(data.frame(ego), 
                    file = paste0(resDir, '/enrichGO_TM3_', cc[n], '.DEgenes.FDR.0.05.txt'),
                    sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
        
        colnames(res1) = paste0(colnames(res1), '_', cc[n], '.pos_vs_RA.pos')
        colnames(res2) = paste0(colnames(res2), '_', cc[n], '.neg_vs_RA.neg')
        res1 = data.frame(res1, res2)
        
      }
      if(n == 1) {
        res = res1
      }else{
        res = data.frame(res, res1)
      }
    }
    
    dev.off()
    
    saveRDS(res, file = paste0(RdataDir, '/TM3_res_pairwiseComparisons.rds'))
    
  }
  
  
  Reduced.model.selection = FALSE
  if(Reduced.model.selection){
    dds1 = dds[, grep('pos', colnames(dds))]
    dds1$conds <- droplevels(dds1$conds)
    dds1 <- estimateDispersions(dds1)
    plotDispEsts(dds1)
    dds1 = nbinomLRT(dds1, reduced = ~ 1)
    res1 <- as.data.frame(results(dds1))
    res1 = res1[order(res1$pvalue), ]
    
    dds2 = dds[, grep('neg', colnames(dds))]
    dds2$conds <- droplevels(dds2$conds)
    dds2 <- estimateDispersions(dds)
    plotDispEsts(dds2)
    dds2 = nbinomLRT(dds2, reduced = ~ 1)
    res2 <- as.data.frame(results(dds2))
    res2 = res2[order(res2$pvalue), ]
    
    
    ggs.signif = unique(c(rownames(res1)[which(res1$pvalue < 0.01 & abs(res1$log2FoldChange) > 0.5)], 
                          rownames(res2)[which(res2$pvalue < 0.001 & abs(res2$log2FoldChange) > 0.5)]))
  }
 
}

##########################################
# edit a list of genes in FGF, BMP and Wnt 
##########################################
Manual.curate.geneList.pathways = FALSE
if(Manual.curate.geneList.pathways){
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

##########################################
# first check the normalized signals in positive and negative cells across conditions  
##########################################
ggs = readRDS(file = paste0(RdataDir, '/TM3_examplesGenes_withGOterm.rds'))
examples = ggs$gene

#examples = ggs.signif
fpm = fpm(dds, robust = TRUE)

#load(file = paste0(RdataDir, '/TM3_positive.negative.pooled_', Counts.to.Use, version.analysis, '.Rdata'))
load(file = paste0(RdataDir, '/TM3_pooled.positive.negative_', Counts.to.Use, version.analysis, '.Rdata'))

fpm.cutoff = 2^2.5

logscale = TRUE
if(logscale){
  fpm = as.matrix(log2(fpm + 2^-4))
  fpm.cutoff = 2.5
  pools = as.matrix(log2(pools + 2^-4))
}

cells = sapply(colnames(fpm), function(x) unlist(strsplit(as.character(x), '_'))[3])
cc = sapply(colnames(fpm), function(x) paste0(unlist(strsplit(as.character(x), '_'))[1:2], collapse = '_'))
cond = sapply(colnames(fpm), function(x) paste0(unlist(strsplit(as.character(x), '_'))[1]))
level_order = apply(expand.grid(c(1:3), c('N2B27', 'RA', 'BMP', 'LDN', 'FGF', 'PD', 'CHIR', 'IWR', 'IWP2')), 
                    1, function(x) paste(x[2], x[1], sep="_"))

#n = which(rownames(fpm) == 'Foxa2')
pdfname = paste0(resDir, '/TM3_examples_FoxA2.positive_negative.pooled_v9_log2scale_pathways.pdf')
pdf(pdfname,  width = 20, height = 16)
#par(cex = 1.0, las = 1, mgp = c(3,2,0), mar = c(6,6,2,0.2), tcl = -0.3)

for(g in examples)
{
  # g = 'Foxa2'
  kk = which(rownames(fpm) == g)
  kk2 = which(rownames(pools) == g)
  
  if(length(kk) > 0){
    cat(g, '\n')
    cpm = fpm[kk, ]
    if(logscale) cpm[which(cpm < 0)] = 0
    xx = data.frame(cpm = cpm, cc = cc,  cells = cells, cond = cond)
    
    p0 = ggplot(xx[which(xx$cells == 'Foxa2.pos'), ],  aes(x = factor(cc, levels = level_order), y = cpm, color = cells, fill = cond)) +
      geom_bar(stat = "identity", position="dodge") +
      geom_hline(yintercept = fpm.cutoff, colour = "blue") + 
      ggtitle(g)  + 
      theme(axis.text.x = element_text(angle = 90, size = 10))
    
    p1 = ggplot(xx[which(xx$cells == 'Foxa2.neg'), ],  aes(x = factor(cc, levels = level_order), y = cpm, color = cells, fill = cond)) +
      geom_bar(stat = "identity", position="dodge") +
      geom_hline(yintercept = fpm.cutoff, colour = "blue") + 
      ggtitle(g)  + 
      theme(axis.text.x = element_text(angle = 90, size = 10))
    
    #grid.arrange(p0, p1, nrow = 1, ncol = 2)
    
    yy = pools[kk2, ]
    if(logscale) yy[which(yy < 0)] = 0
    yy = data.frame(cpm = yy, cc = colnames(pools), cells = rep('pooled', length(yy)), cond = cc.pools)
    p2 = ggplot(yy,  aes(x = factor(cc, levels = level_order), y = cpm, color = cells, fill = cond)) +
      geom_bar(stat = "identity", position="dodge") +
      geom_hline(yintercept = fpm.cutoff, colour = "blue") +
      ggtitle(g)  +
      theme(axis.text.x = element_text(angle = 90, size = 10))
    
    grid.arrange(p0, p1, p2, nrow = 2, ncol = 2)
    
  }
  
}

dev.off()


########################################################
########################################################
# Section : test pooling of positive and negative cells
# 
########################################################
########################################################
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

##########################################
# explore other players 
##########################################

