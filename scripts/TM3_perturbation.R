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
RdataDir = paste0(resDir, '/Rdata/')

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

##################################################
## Import UMI count table
##################################################
source(RNA.functions)
source(RNA.QC.functions)

Import.onlyUMI = TRUE

design = readRDS(file = paste0(RdataDir, 'sampleInfo_QC.stats.rds'))

colnames(design)[1] = 'SampleID'
design$conds = paste0(design$condition, '_', design$cells)

# prepare the data table 
Dir_umi = paste0(dataDir, "htseq_counts_BAMs_umi")
Dir_read = paste0(dataDir, "htseq_counts_BAMs")

aa1 <- list.files(path = Dir_umi, pattern = "*umiDedup.txt", full.names = TRUE)
aa1 = merge.countTables.htseq(aa1)
colnames(aa1)[-1] = paste0(colnames(aa1)[-1], ".UMI")

if(!Import.onlyUMI){
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
  
}else{
  aa = aa1
}

# if(Counts.to.Use == 'readCounts'){
#    all = process.countTable(all=aa, design = design, special.column = ".readCount", ensToGeneSymbol = FALSE)
#  }else{
#    if(Counts.to.Use == "UMI"){
#      all = process.countTable(all=aa, design = design, special.column = "UMI", ensToGeneSymbol = FALSE)
#    }else{
#      cat("Error : no counts found for ", Counts.to.Use, "for miRNAs \n")
#    }
#  }


all = process.countTable(all=aa, design = design[,c(1, 14)], special.column = NULL, ensToGeneSymbol = FALSE)

all = all[grep('^__', all$gene, invert = TRUE), ]

save(design, all, file=paste0(RdataDir, 'Design_Raw_readCounts', Counts.to.Use,  version.analysis, '.Rdata'))

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
    
    mm = match(all$gene, annot$Gene.stable.ID)
    cat(length(which(!is.na(mm))), ' protein coding genes in the count table \n')
    
    all = all[!is.na(mm), ]
    
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

save(design, all, file = paste0(RdataDir, 'design.detailed_RawUMI_', Counts.to.Use, version.analysis, '.Rdata'))

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
# Section : normalization with DESeq2 
# 
########################################################
########################################################
require(DESeq2)
require(ggplot2)
library(ggrepel)
require(gridExtra)
library(dplyr)
library(patchwork)
require(pheatmap)
library(org.Mm.eg.db)
library(enrichplot)
library(clusterProfiler)
library(stringr)

Counts.to.Use = 'UMI'

load(file = paste0(RdataDir, '/design.detailed_RawUMI_', Counts.to.Use, version.analysis, '.Rdata'))
genes = readRDS(file = paste0(RdataDir, '/mm10_ens_BioMart_GRCm38.p6_ensID.geneSymbol.length.rds'))
ggs = readRDS(file = paste0('../results/Rdata/curated_signaling.pathways_gene.list_v2.rds'))

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

length(which(ss > quantile(ss, probs = 0.6)))
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

save(dds, design.matrix,  file = paste0(RdataDir, 'TM3_dds_normalized.Rdata'))
#save(fpm, design, file = paste0(tfDir, '/RNAseq_fpm_fitered.cutoff.', cutoff.gene, '.Rdata'))
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

ss = apply(fpm, 1, mean)
jj = which(ss < quantile(ss, probs = 0.3))
length(jj)

pca=plotPCA(vsd[jj, ], intgroup = c('condition', 'cells'), returnData = TRUE, ntop = 500)
pca2save = as.data.frame(pca)

ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition, shape = cells))  + 
  geom_point(size=4) + 
  geom_text(hjust = 0.2, nudge_y = 0.25, size=3)

plot(ggp) + ggsave(paste0(resDir, "/PCAplot_withControls_UQ.norm_ntop500.pdf"), width=18, height = 12)

ggp = ggplot(data=pca2save[grep('N2B27', pca2save$condition, invert = TRUE), ], 
             aes(PC1, PC2, label = name, color= condition, shape = cells))  + 
  geom_point(size=4) + 
  geom_text(hjust = 0.2, nudge_y = 0.1, size=3)

plot(ggp) + ggsave(paste0(resDir, "/PCAplot_withoutControls_UQ.norm_ntop500.pdf"), width=18, height = 12)

##########################################
# try to put the TM3 data in the context of time series data
# test the relationship between perturbation of positive and negative cells, pooled cells, and time seires and day 5 sorted cells 
# 
##########################################
Compare.with.pooledCells.timeSeries.sortedDay5 = FALSE
if(Compare.with.pooledCells.timeSeries.sortedDay5){
  source('Functions_RNAseqData.R')
  
  Test.pooling.TM3.positive.negative.cells(dds, pooling = 'counts', design.matrix = design.matrix)
  
  Compare.TM3.and.RNAseq.timeSeries.sortedDay5(dds)
  
  Compare.TM3.sortedDay3.RNAseq.sortedDay5(dds)
  
}

########################################################
########################################################
# Section : Pairwise comparisons 
# Here we first identify significantly different genes in RA positive and negative genes
# RA positive vs negative
# In addition, the different response of Foxa2 positive and negative cells are also containing information of 
# which genes are expressed in positive or negative cells
########################################################
########################################################
Compare.RA.positve.negative = FALSE
if(Compare.RA.positve.negative){
  
  kk = which(design.matrix$condition == 'RA')
  
  dds1 = dds[, kk]
  
  dds1$conds <- droplevels(dds1$conds)
  
  #dds1 <- estimateDispersions(dds1, fitType = 'parametric')
  #plotDispEsts(dds1, ymin = 10^-3, main = 'RA')
  sfs = sizeFactors(dds1)
  dds1 = DESeqDataSetFromMatrix(counts(dds1), DataFrame(design.matrix[kk, ]), design = ~ conds + condition.replicate)
  sizeFactors(dds1) = sfs
  
  dds1 <- estimateDispersions(dds1, fitType = 'parametric')
  plotDispEsts(dds1, ymin = 10^-3, main = 'RA')
  
  dds1 <- nbinomWaldTest(dds1)
  resultsNames(dds1)  
  
  res1 = results(dds1, contrast=c("conds", 'RA_Foxa2.pos', 'RA_Foxa2.neg'), alpha = 0.05)
  res1 <- lfcShrink(dds1, coef=2, type="normal")
  
  plotMA(res1)
  
  res1 = as.data.frame(res1)
  gg.signif = rownames(res1)[which(res1$padj < 0.05)]
  
  # save the average expression of Foxa2 positive and negative cells
  cpm = fpm(dds1)
  cpm = cpm[, c(grep('pos', colnames(cpm)), grep('neg', colnames(cpm)))]
  
  res1 = data.frame(res1, cpm.foxa2.pos = apply(as.matrix(cpm[, grep('pos', colnames(cpm))]), 1, mean), 
                    cpm.foxa2.neg = apply(as.matrix(cpm[, grep('neg', colnames(cpm))]), 1, mean))
  
  saveRDS(res1, file = paste0(resDir, '/TM3_Foxa2.positive_vs_neg.day3_comparison_cpm.mean.rds'))
  
  save(res1, cpm, file = paste0(resDir, '/TM3_Foxa2.positive_vs_neg.day3_comparisonStats_cpm.Rdata'))
  
  ##########################################
  # ## visualize global DE genes 
  #  ## GO enrich analysis for upregulated and downregulated genes
  ##########################################
  df <- data.frame(cond = colnames(cpm), sorted = c(rep('Foxa2.pos', 3), rep('Foxa2.neg', 3)))
  rownames(df) = colnames(cpm)
  
  pheatmap(cpm[which(res1$padj<0.05), ], cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = TRUE, 
           cluster_cols=FALSE, na_col = 'gray',  annotation_col = df, gaps_col = c(3),
           filename = paste0(resDir, '/TM3_RAtreatment_DEgenes_Day3.pdf'), width = 10, height = 8)
  
  xx = data.frame(res1[order(res1$pvalue), ])
  
 
  bgs = rownames(dds1)
  bgs.df <- bitr(bgs, fromType = "SYMBOL",
                 toType = c("ENSEMBL", "ENTREZID"),
                 OrgDb = org.Mm.eg.db)
  
  gg.signif = rownames(res1)[which(res1$padj < 0.1 & res1$log2FoldChange<0)]
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
  
  p0 = barplot(ego, showCategory=50, main = 'RA') + 
    ggsave(paste0(resDir, "/TM3_RA.treatment_positive.vs.negative.DEgenes.downregulated_top50.pdf"), width=12, height = 12)
  ii = grep('pathway', ego[, 2])
  xx = ego[ii, ]
  par(cex = 1.0, las = 1, mgp = c(3,1,0), mar = c(3, 30,2,0.2), tcl = -0.3)
  barplot(xx$Count,  names.arg = xx$Description, horiz = TRUE, beside = TRUE, las=2)
  
  xx = data.frame(res1[order(res1$pvalue), ])
  xx = xx[which(xx$padj < 0.05), ]
  xx = xx[!is.na(match(rownames(xx), ggs$gene)), ]
  
  ##########################################
  # plot the fold-changes for Wnt, FGF, BMP siganling pathways
  ##########################################
  library(tidyr)
  library(dplyr)
  res = readRDS(file = paste0(resDir, '/TM3_Foxa2.positive_vs_neg.day3_comparison_cpm.mean.rds'))
  load(file = paste0('../results/Rdata/RNAseq_Foxa.positive_vs_neg.Day5.Rdata'))
  gene.all = unique(intersect(rownames(res), rownames(res1)))
  
  res = data.frame(res[match(gene.all, rownames(res)), c(2, 5, 7, 8)], res1[match(gene.all, rownames(res1)), c(2, 5, 7, 8)])
  colnames(res) = c('lfc.day3', 'pval.day3', 'cpm.pos.day3', 'cpm.neg.day3',
                    'lfc.day5', 'pval.day5', 'cpm.pos.day5', 'cpm.neg.day5')
  
  res$gene = rownames(res)
  
  #yy1 = cbind(res[,c(1, 2, 5)], rep('day3', nrow(res)))
  #yy2 = cbind(res[,c(3, 4, 5)], rep('day5', nrow(res)))
  #colnames(res1) = paste0(colnames(res1), '_RA.pos_vs_RA.neg')
  #res = readRDS(file = paste0(RdataDir, '/TM3_res_pairwiseComparisons.rds'))
  #examples = c('Foxa2', 'Olig2', 'Shh')
  
  ##########################################
  # WNT, FGF, BMP genes in Foxa2 positive and negative cells
  ##########################################
  examples = c('Lef1', 'Wnt3', 'Wnt3a', 'Wnt4', 'Wnt5b', 'Wnt6', 'Wnt7a', 'Wnt7b', 'Wnt8a', 'Dkk1', 'Dkk2', 'Dkk3', 'Tcf15', 'Tcf19', 
               "Wnt1", 'Sost', 'Sfrp5', 'Lypd6', 'Peg12', 'Notum', 'Draxin')
  SP = 'WNT'
  
  examples = c('Spry4', 'Spry2', 'Etv4', 'Etv5', 'Fgf10', 'Fgf17','Fgf8', 'Fgf5', 'Fgf2','Fgf21', 'Fgf11','Fgf1', 'Fgf4', 'Fgfbp3',
               'Dusp1', 'Dusp10', 'Dusp27', 'Dusp4', 'Dusp5')
  SP = 'FGF'
  
  examples = c('Id1', 'Id3', 'Smad6', 'Nog', 'Fst', 'Bambi', 'Bmp7', 'Bmp4', 'Bmp1', 'Bmp6', 'Bmpr2', 'Bmpr1b', 'Bmp2', 'Bmp3', 'Bmp5', 
               'Runx1', 'Runx2')
  SP = 'BMP'
  
    
  kk = which(!is.na(match(rownames(res), examples)) == TRUE)
  yy = res[, c(1, 3, 4, 9)]; colnames(yy) = c('lfc', 'cpm.pos', 'cpm.neg', 'gene')
  yy$regulation = 'Foxa2+'; yy$regulation[which(yy$lfc<0)] = 'Foxa2-'
  
  ggplot(yy[kk, ], aes(x = gene, y = lfc, fill = regulation)) + 
    geom_bar(stat = "identity", position="dodge") + ylim(-1, 1) +
    coord_flip() + 
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
    labs(title = paste0(SP, ' day3'), y = 'LFC', x = '') + 
    ggsave(paste0(resDir, "/LFC_positive.vs.negative_Day3_", SP, ".pdf"), width=8, height = 6)
  
  as_tibble(yy[kk, ]) %>% gather(cells, cpm, 2:3) %>%
    ggplot(aes(x = gene, y = cpm, fill = cells)) +
    geom_bar(stat = "identity", position="dodge", width = 0.7) + 
    coord_flip() +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
    labs(title = paste0(SP, ' day3'), y = 'CPM', x = '') + 
    ggsave(paste0(resDir, "/CPM_positive.vs.negative_Day3_", SP, ".pdf"), width=8, height = 10)
  
  
  yy = res[, c(5, 7, 8, 9)]; colnames(yy) = c('lfc', 'cpm.pos', 'cpm.neg', 'gene')
  yy$regulation = 'Foxa2+'; yy$regulation[which(yy$lfc<0)] = 'Foxa2-'
  
  ggplot(yy[kk, ], aes(x = gene, y = lfc, fill = regulation)) + 
    geom_bar(stat = "identity", position="dodge")  +
    coord_flip() + 
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
    labs(title = paste0(SP, ' day5'), y = 'LFC', x = '') + 
    ggsave(paste0(resDir, "/LFC_positive.vs.negative_Day5_", SP, ".pdf"), width=8, height = 6)
  
  as_tibble(yy[kk, ]) %>% gather(cells, cpm, 2:3) %>%
    ggplot(aes(x = gene, y = cpm, fill = cells)) +
    geom_bar(stat = "identity", position="dodge", width = 0.7) + 
    coord_flip() +
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12)) + 
    labs(title = paste0(SP,  ' day5'), y = 'CPM', x = '') + 
    ggsave(paste0(resDir, "/CPM_positive.vs.negative_Day5_", SP, ".pdf"), width=8, height = 10)
  
}

########################################################
########################################################
# Section : Pairwise comparisons to identify potential interactions 
# Here we compare 
# -perturbed-positive vs RA-positive
# -perturbed-negative vs RA-negative
# -perturbed-pooled vs RA-pooled
########################################################
########################################################
Calculate.pairwise.comparisons = FALSE
Compare.pooled.samples = TRUE

if(Calculate.pairwise.comparisons){
  if(Compare.pooled.samples){
    # save(ddx, cc.pools, file = paste0(RdataDir, '/TM3_pooled.pos.neg_ddx_cc.pools_', Counts.to.Use, version.analysis, '.Rdata'))
    load(file = paste0(RdataDir, '/TM3_pooled.pos.neg_ddx_cc.pools_', Counts.to.Use, version.analysis, '.Rdata'))
    
  }
  
  cc = unique(design.matrix$condition)
  cc = cc[which(cc != "N2B27" & cc != 'RA')]
  
  bgs = rownames(dds)
  bgs.df <- bitr(bgs, fromType = "SYMBOL",
                 toType = c("ENSEMBL", "ENTREZID"),
                 OrgDb = org.Mm.eg.db)
  
  pdfname = paste0(resDir, '/TM3_pairwiseComparisons_enrichGO_pathways_sorted.positve.negative.pooled.pdf')
  pdf(pdfname,  width = 18, height = 16)
  #par(cex = 1.0, las = 1, mgp = c(3,2,0), mar = c(6,6,2,0.2), tcl = -0.3)
  
  for(n in 1:length(cc))
  {
    # n = 1
    cat(as.character(cc[n]), '\n')
    
    kk = which(design.matrix$condition == cc[n] | design.matrix$condition == 'RA')
    dds1 = dds[, kk]
    dds1$conds <- droplevels(dds1$conds)
    dds1 <- estimateDispersions(dds1, fitType = 'parametric')
    plotDispEsts(dds1, ymin = 10^-3, main = cc[n])
    #dds = estimateDispersions(dds)
    #plotDispEsts(dds)
    dds1 <- nbinomWaldTest(dds1)
    resultsNames(dds1)  
    
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
    
    ## add the comparison of pooled samples
    if(Compare.pooled.samples){
      dds2 = ddx[, which(cc.pools == cc[n] | cc.pools == 'RA')]
      dds2$condition <- droplevels(dds2$condition)
      dds2 <- estimateDispersions(dds2, fitType = 'parametric')
      plotDispEsts(dds1, ymin = 10^-3, main = paste0(cc[n], ' pooled'))
      #dds = estimateDispersions(dds)
      #plotDispEsts(dds)
      dds2 <- nbinomWaldTest(dds2)
      resultsNames(dds2)  
      
      res3 = results(dds2, contrast=c("condition", cc[n], 'RA'), alpha = 0.05)
      res3 <- lfcShrink(dds2, contrast=c("condition", cc[n], 'RA'), type="normal")
      res3 = as.data.frame(res3)
      
      colnames(res3) = paste0(colnames(res3), '_', cc[n], '.pooled_vs_RA.pooled')
      
      res1 = data.frame(res1, res3)
      
    }
    
    ii = grep('log2FoldChange', colnames(res1))
    plot.pair.comparison.plot(res1[, ii], linear.scale = FALSE)
    
    ii = grep('pvalue', colnames(res1))
    plot.pair.comparison.plot(-log10(res1[, ii]), linear.scale = FALSE)
    
    if(n == 1) {
      res = res1
    }else{
      res = data.frame(res, res1)
    }
  }
  
  dev.off()
  
  saveRDS(res, file = paste0(RdataDir, '/TM3_res_pairwiseComparisons_perturbation.vs.RA_positive.negative.pooled.rds'))
  
  ##########################################
  # extract enriched signaling pathways from perturbations
  ##########################################
  Extract.signaling.pathways.enrichedGOterm = FALSE
  if(Extract.signaling.pathways.enrichedGOterm){
    files = list.files(path = resDir, pattern = 'DEgenes.FDR.0.05.txt', full.names = TRUE)
    files = files[c(7, 1:6)]
    
    for(n in 1:length(files))
    {
      # n = 4
      ff = read.table(files[n], sep = '\t', header = TRUE, row.names = 1, as.is = c(1, 2,8))
      
      #ids = rbind(ids, cbind(ff$ID, as.character(ff$Description), as.character(ff$geneID) ))
      ff = data.frame(ff, stringsAsFactors = FALSE)
      #colnames(ids) = c('ID', 'description', 'geneID')
      ii = grep('pathway', ff[, 2])
      xx = ff[ii, ]
      xx = xx[order(-xx$p.adjust), ]
      xx$Description = gsub('signaling', '', xx$Description)
      xx$Description = gsub(' pathway', '', xx$Description)
      par(cex = 1.0, las = 1, mgp = c(3,1,0), mar = c(3, 30,2,0.2), tcl = -0.3)
      barplot(xx$Count,  names.arg = xx$Description, horiz = TRUE, beside = TRUE, las=2)
      
    }
    
  }
  
  ##########################################
  # Here each perturbation will be checked for Foxa2 positive, negative and pooled
  # Top response genes will be highlighted, some of FGF, Wnt and BMP related genes will be highlighted as well
  # the p-values and lfc will be compared between positve, negative and pooled cells
  ##########################################
  Explore.perturbation.response.for.signaling.pathway.genes = FALSE
  if(Explore.perturbation.response.for.signaling.pathway.genes){
    
    res = readRDS(file = paste0(RdataDir, '/TM3_res_pairwiseComparisons_perturbation.vs.RA_positive.negative.pooled.rds'))
    ggs = readRDS(file = paste0(RdataDir, '/TM3_examplesGenes_withGOterm.rds'))
    
    Compare.positive.vs.pooled.samples = FALSE
    
    cc = unique(design.matrix$condition)
    cc = cc[which(cc != "N2B27" & cc != 'RA')]
    
    res = res[!is.na(match(rownames(res), ggs$gene)), grep('pvalue_|log2FoldChange_', colnames(res))]
    
    library(ggrepel)
    library(ggplot2)
    
    examples = c('Bmp4', 'Bmp1', 'Bmp7', 'Bmp6', 'Nog', 'Dkk1', 'Fgfbp3', 'Lef1', 'Tcf15', 'Tcf19', 'Wnt3', 'Wnt3a', 'Wnt4', 
                 'Wnt5b', 'Wnt6', 'Wnt7a', 'Wnt7b', 'Wnt8a', 'Sfrp5', 'Lypd6', 'Spry4', 'Etv5', 'Etv4', 'Fgf10', 'Fgf17', 
                 'Fgf8', 'Fgf5', 'Fgf2', 'Fgf21', 'Fgf11', 'Fgf1', 'Fgf4')
    
    
    
    pdfname = paste0(resDir, '/TM3_comparing_response_signalingPathways_positve.negative_v2.pdf')
    pdf(pdfname,  width = 16, height = 12)
    #par(cex = 1.0, las = 1, mgp = c(3,2,0), mar = c(6,6,2,0.2), tcl = -0.3)
    
    for(n in 1:length(cc))
    {
      # n = 5
      cat(n, ' -- ', cc[n], '\n')
      jj = grep(cc[n], colnames(res))
      yy = res[, jj]
      colnames(yy) = gsub(paste0(cc[n]), '', colnames(yy))
      colnames(yy) = gsub('_vs_RA.', '', colnames(yy))
      colnames(yy) = gsub('.pos|.neg|.pooled', '', colnames(yy))
      colnames(yy) = gsub('log2FoldChange', 'lfc', colnames(yy))
      
      yy[, grep('pvalue_', colnames(yy))] = -log10(yy[, grep('pvalue_', colnames(yy))])
      yy = data.frame(yy, gene = rownames(yy))
      
      p1 = ggplot(data = yy,  aes(y = pvalue_pos, x = lfc_pos,  label = gene)) + 
        geom_point(size = 1) + 
        labs(title = paste0(cc[n], " - positive  "), x = '', y = '-log10(pval)') + 
        theme(axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12)) + 
        geom_hline(yintercept = 2, colour = "red") +
        #geom_text(data=subset(yy, pvalue_pos > 2), size = 4, nudge_y = 0.5) + 
        geom_text_repel(data=subset(yy, pvalue_pos > 2), size = 4)
      
      p2 = ggplot(data = yy,  aes(y = pvalue_neg, x = lfc_neg,  label = gene)) + 
        geom_point(size = 1) + 
        labs(title = paste0(cc[n], " - neg  "), x = '', y = '-log10(pval)') + 
        theme(axis.text.x = element_text(size = 12), 
              axis.text.y = element_text(size = 12)) + 
        geom_hline(yintercept = 2, colour = "red") +
        #geom_text(data=subset(yy, pvalue_pos > 2), size = 4, nudge_y = 0.5) + 
        geom_text_repel(data=subset(yy, pvalue_neg > 2), size = 4)
      
      if(Compare.positive.vs.pooled.samples){
        p3 = ggplot(data = yy,  aes(y = pvalue_pooled, x = lfc_pooled,  label = gene)) + 
          geom_point(size = 1) + 
          labs(title = paste0(cc[n], " - pooled  "), x = '', y = '-log10(pval)') + 
          theme(axis.text.x = element_text(size = 12), 
                axis.text.y = element_text(size = 12)) + 
          geom_hline(yintercept = 2, colour = "red") +
          #geom_text(data=subset(yy, pvalue_pos > 2), size = 4, nudge_y = 0.5) + 
          geom_text_repel(data=subset(yy, pvalue_pooled > 2), size = 4)
      }
     
      examples.sel = unique(c(examples, rownames(yy)[which(yy$pvalue_pos > 2 | yy$pvalue_neg > 2)]))
      
      if(Compare.positive.vs.pooled.samples){
        p4 = ggplot(data = yy, aes(x = pvalue_pooled, y = pvalue_neg , label = gene)) +
          geom_point(size = 1) + 
          xlim(0, 6) + ylim(0, 6) + geom_hline(yintercept = 1.3, colour = "blue") +
          geom_vline(xintercept = 1.3, colour = "blue") + 
          labs(title = paste0(cc[n], " - neg vs pooled : -log10(pval) "), x = 'pooled', y = 'negative') +
          geom_abline(slope = 1, intercept = 0, colour = 'red') + 
          geom_label_repel(data=  as.tibble(yy) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel),
                           size = 3)
        #examples.sel = unique(c(examples, rownames(yy)[which(yy$pvalue_pos > 2 | yy$pvalue_neg > 2)]))
        p5 = ggplot(data = yy, aes(x = pvalue_pooled, y = pvalue_pos, label = gene)) +
          geom_point(size = 1) + 
          xlim(0, 6) + ylim(0, 6) + geom_hline(yintercept = 1.3, colour = "blue") +
          geom_vline(xintercept = 1.3, colour = "blue") + 
          labs(title = paste0(cc[n], " - pos vs pooled : -log10(pval) "), x = 'pooled', y = 'positive') + 
          geom_abline(slope = 1, intercept = 0, colour = 'red') + 
          geom_label_repel(data=  as.tibble(yy) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel),
                           size = 3)
      }
      
      p6 = ggplot(data = yy, aes(x = pvalue_neg, y = pvalue_pos, label = gene)) +
        geom_point(size = 1) + 
        xlim(0, 6) + ylim(0, 6) + geom_hline(yintercept = 2, colour = "blue") +
        geom_vline(xintercept = 2, colour = "blue") + 
        labs(title = paste0(cc[n], " - pos vs neg : -log10(pval) "), x = 'negative', y = 'positive') + 
        geom_abline(slope = 1, intercept = 0, colour = 'red') + 
        geom_label_repel(data=  as.tibble(yy) %>%  dplyr::mutate_if(is.factor, as.character) %>% dplyr::filter(gene %in% examples.sel),
                         size = 3)
      
      plot(p1)
      plot(p2)
      if(Compare.positive.vs.pooled.samples){
        plot(p3)
        plot(p4)
        plot(p5)
      }
      plot(p6)
      
      #grid.arrange(p1, p2, p3, p4, p5, nrow = 2, ncol = 3)
    }
    
    dev.off()
    
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


########################################################
########################################################
# Section : plotting results
# 1) pairwise comaprison between RA positive and RA negative
# 2) # first check the normalized signals in positive and negative cells across conditions  
########################################################
########################################################
##########################################
# plot individual examples
##########################################
load(file = paste0(RdataDir, 'TM3_dds_normalized.Rdata'))
#ggs = readRDS(file = paste0(RdataDir, '/TM3_examplesGenes_withGOterm.rds'))
ggs = readRDS(file = paste0('../results/Rdata/curated_signaling.pathways_gene.list_v2.rds'))

ggs$pathway[which(is.na(ggs$pathway))] = 'A' 

#ggs = ggs[order(ggs$pathway), ]
ggs = ggs[which(ggs$pathway == 'SHH'), ]

res = readRDS(file = paste0(RdataDir, '/TM3_res_pairwiseComparisons_perturbation.vs.RA_positive.negative.pooled.rds'))
#load(file = paste0(RdataDir, '/TM3_positive.negative.pooled_', Counts.to.Use, version.analysis, '.Rdata'))
#load(file = paste0(RdataDir, '/TM3_pooled.positive.negative_', Counts.to.Use, version.analysis, '.Rdata'))
load(file = paste0(RdataDir, '/TM3_pooled.pos.neg_ddx_cc.pools_', Counts.to.Use, version.analysis, '.Rdata'))

examples = ggs$gene
pvals = res[, grep('pvalue_', colnames(res))]
p.min = apply(as.matrix(pvals), 1, min)
kk = which(p.min<0.05)
examples = examples[!is.na(match(examples, rownames(res)[kk]))]

#examples = ggs.signif
fpm = fpm(dds, robust = TRUE)
pools = fpm(ddx, robust = TRUE)

fpm = fpm[, grep('N2B27', colnames(fpm), invert = TRUE)]
pools = pools[, grep('N2B27', colnames(pools), invert = TRUE)]

fpm.cutoff = 2^2.5

logscale = TRUE
if(logscale){
  fpm = as.matrix(log2(fpm + 2^-4))
  fpm.cutoff = 2.
  pools = as.matrix(log2(pools + 2^-4))
}


cells = sapply(colnames(fpm), function(x) unlist(strsplit(as.character(x), '_'))[3])
cc = sapply(colnames(fpm), function(x) paste0(unlist(strsplit(as.character(x), '_'))[1:2], collapse = '_'))
cond = sapply(colnames(fpm), function(x) paste0(unlist(strsplit(as.character(x), '_'))[1]))
level_order = apply(expand.grid(c(1:3), c('RA', 'BMP', 'LDN', 'FGF', 'PD', 'CHIR', 'IWP2', 'IWR')), 
                    1, function(x) paste(x[2], x[1], sep="_"))

#n = which(rownames(fpm) == 'Foxa2')
pdfname = paste0(resDir, '/TM3_examples_FoxA2.positive_negative.pooled_v12_log2scale_signalingPathwaySHH.pdf')
pdf(pdfname,  width = 20, height = 16)
#par(cex = 1.0, las = 1, mgp = c(3,2,0), mar = c(6,6,2,0.2), tcl = -0.3)
library(reshape2)
for(g in examples)
{
  # g = 'Id1'
  kk = which(rownames(fpm) == g)
  kk2 = which(rownames(pools) == g)
  kk3 = which(rownames(res) == g)
  if(length(kk) > 0){
    cat(g, '\n')
    cpm = fpm[kk, ]
    if(logscale) cpm[which(cpm < -2)] = -2
    xx = data.frame(cpm = cpm, cc = cc,  cells = cells, cond = cond)
    
    p0 = ggplot(xx[which(xx$cells == 'Foxa2.pos'), ],  aes(x = factor(cc, levels = level_order), y = cpm, fill = cond)) +
      geom_bar(stat = "identity", position="dodge") +
      geom_hline(yintercept = fpm.cutoff, colour = "blue") + 
      labs(title = paste0(g, " - positive  "), x = '', y = 'log2(cpm)') + 
      theme(axis.text.x = element_text(angle = 90, size = 10))
    
    p1 = ggplot(xx[which(xx$cells == 'Foxa2.neg'), ],  aes(x = factor(cc, levels = level_order), y = cpm, fill = cond)) +
      geom_bar(stat = "identity", position="dodge") +
      geom_hline(yintercept = fpm.cutoff, colour = "blue") + 
      labs(title = paste0(g, " - negative "), x = '', y = 'log2(cpm)') + 
      theme(axis.text.x = element_text(angle = 90, size = 10))
    
    #grid.arrange(p0, p1, nrow = 1, ncol = 2)
    yy = pools[kk2, ]
    if(logscale) yy[which(yy < -2)] = -2
    yy = data.frame(cpm = yy, cc = colnames(pools), cells = rep('pooled', length(yy)), cond = cc.pools[which(cc.pools != 'N2B27')])
    p2 = ggplot(yy,  aes(x = factor(cc, levels = level_order), y = cpm,  fill = cond)) +
      geom_bar(stat = "identity", position="dodge") +
      geom_hline(yintercept = fpm.cutoff, colour = "blue") +
      labs(title = paste0(g, " - pooled  "), x = '', y = 'log2(cpm)') + 
      theme(axis.text.x = element_text(angle = 90, size = 10))
    
    zz = res[kk3, grep('pvalue_', colnames(res))]
    zz = matrix(as.numeric(zz), nrow = 3, ncol = 7)
    rownames(zz) = c('positive', 'negative', 'pooled')
    cz = unique(design.matrix$condition)
    cz = cz[which(cz != "N2B27" & cz != 'RA')]
    colnames(zz) = cz
    zz = -log10(zz)
    
    zz[which(zz>6)] = 6
    
    zz_melted<-melt(zz)
    
    p3 = ggplot(zz_melted, aes(x = Var2, y = factor(Var1, levels = c('pooled', 'negative',  'positive')))) +  
      geom_point( aes(size=value),shape=21, colour="black", fill="blue") +
      theme(panel.background=element_blank(), panel.border = element_rect(colour = "blue", fill=NA, size=2), 
            axis.text=element_text(size=12)) +
      scale_size_area(max_size = 12, breaks = c(2, 4, 6)) +
      #Add labels to axes
      labs(x="", y="", title = '-log10(p-val)')
    
    grid.arrange(p0, p1, p2, p3, nrow = 2, ncol = 2)
    
  }
  
}

dev.off()


##########################################
# make plots for Noggin in BMP perturbation 
##########################################
Make.plot.examples.noggin = FALSE
if(Make.plot.examples.noggin){
  
  load(file = paste0(resDir, '/TM3_Foxa2.positive_vs_neg.day3_comparisonStats_cpm.Rdata'))
  
  res = readRDS(file = paste0(RdataDir, '/TM3_res_pairwiseComparisons_perturbation.vs.RA_positive.negative.pooled.rds'))
  
  load(file = paste0(RdataDir, '/TM3_pooled.pos.neg_ddx_cc.pools_', Counts.to.Use, version.analysis, '.Rdata'))
  fpm = fpm(dds, robust = TRUE)
  gg = 'Nog'
  fpm = fpm[which(rownames(fpm) == gg), grep('RA_|BMP_', colnames(fpm))]
  fpm = log2(fpm)
  fpm = data.frame(data = fpm)
  fpm$condition = sapply(rownames(fpm), function(x) unlist(strsplit(as.character(x), '_'))[1])
  fpm$cell = sapply(rownames(fpm), function(x) unlist(strsplit(as.character(x), '_'))[3])
  
  fpm$cc = paste0(fpm$condition, '_', fpm$cell)
  fpm$sd = NA
  fpm$mean = NA
  for(n in 1:nrow(fpm))
  {
    jj = which(fpm$cc == fpm$cc[n])
    fpm$mean[n] = mean(fpm$data[jj])
    fpm$sd[n] = sd(fpm$data[jj])/sqrt(length(jj))
  }
  
  fpm = fpm[match(unique(fpm$cc), fpm$cc), ]
  
  library(ggplot2)
  # Default bar plot
  p<- ggplot(fpm, aes(x=factor(condition, levels = c('RA', 'BMP')), y=mean, fill= factor(cell, levels = c('Foxa2.pos', 'Foxa2.neg')))) + 
    geom_bar(stat="identity", color="black", 
             position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(.9)) +
    theme(legend.title=element_blank()) +
    labs(x="", y="log2(normalized data) ")
  print(p) + ggsave(paste0(resDir, "/Noggin_example.pdf"), width=8, height = 8)
  
  
  
}


