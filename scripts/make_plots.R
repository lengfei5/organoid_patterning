##########################################################################
##########################################################################
# Project: Teresa's organoid project
# Script purpose: make plots for the paper
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Tue Jun 13 17:20:43 2023
##########################################################################
##########################################################################
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
# import the data and results
##########################################
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

load(file = paste0(RdataDir, 'TM3_dds_normalized.Rdata'))

#load(file = paste0(RdataDir, '/design.detailed_RawUMI_', Counts.to.Use, version.analysis, '.Rdata'))
#load(file = paste0(RdataDir, '/TM3_pooled.pos.neg_ddx_cc.pools_', Counts.to.Use, version.analysis, '.Rdata'))

cc = unique(design.matrix$condition)
cc = c('RA', 'BMP', 'LDN')


kk = which(!is.na(match(design.matrix$condition, cc)))
design = design.matrix[kk, ]
dds1 = dds[, kk]
dds1$conds <- droplevels(dds1$conds)
dds1 <- estimateDispersions(dds1, fitType = 'parametric')
plotDispEsts(dds1, ymin = 10^-3)

dds1 <- nbinomWaldTest(dds1)
resultsNames(dds1)  

cpm = fpm(dds1, robust = TRUE)

compares = c()
for(n in seq_len(length(cc)))
{
  compares = c(compares, list(c('conds', paste0(cc[n], '_Foxa2.pos'), paste0(cc[n], '_Foxa2.neg'))))
}

compares = c(compares, 
             list(c('conds', paste0('BMP', '_Foxa2.pos'), paste0('RA', '_Foxa2.pos'))), 
             list(c('conds', paste0('BMP', '_Foxa2.neg'), paste0('RA', '_Foxa2.neg'))),
             list(c('conds', paste0('LDN', '_Foxa2.pos'), paste0('RA', '_Foxa2.pos'))),
             list(c('conds', paste0('LDN', '_Foxa2.neg'), paste0('RA', '_Foxa2.neg'))),
             list(c('conds', paste0('BMP', '_Foxa2.pos'), paste0('LDN', '_Foxa2.pos'))),
             list(c('conds', paste0('BMP', '_Foxa2.neg'), paste0('LDN', '_Foxa2.neg'))))

ggs = c()
outDir = paste0(resDir, '/plots/')
if(!dir.exists(outDir)) dir.create(outDir)

for(n in seq_len(length(compares)))
{
  # n = 4
  res = results(dds1, contrast=compares[[n]], alpha = 0.1)
  res<- lfcShrink(dds1, contrast=compares[[n]], type="normal")
  res = as.data.frame(res)
  
  # save the significant genes
  gg.signif = rownames(res)[which(res$padj < 0.1 & abs(res$log2FoldChange) > 1)] 
  cat(length(gg.signif), ' significant genes for comparison: ', compares[[n]], '  \n')
  ggs = c(ggs, gg.signif)
  
  yy = res
  colnames(yy) = gsub('log2FoldChange', 'lfc', colnames(yy))
  yy$pvalue = -log10(yy$padj)
  yy = data.frame(yy, gene = rownames(yy))
  yy = yy[which(!is.na(yy$lfc & !is.na(yy$pvalue))), ]
  
  ggplot(data = yy,  aes(y = pvalue, x = lfc,  label = gene)) + 
    geom_point(size = 1) + 
    labs(title = paste0(compares[[n]][c(2:3)], collapse = '_vs._'), x = 'log2FC', y = '-log10(pval)') + 
    theme_classic()+
    theme(axis.text.x = element_text(size = 12), 
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12)) + 
    #geom_hline(yintercept = 2, colour = "red") +
    #geom_text(data=subset(yy, pvalue_pos > 2), size = 4, nudge_y = 0.5) + 
    geom_text_repel(data=subset(yy, pvalue > 3|gene == 'Nog'), size = 4) 
  
  ggsave(filename = paste0(outDir, 'volcanoPlot_', paste0(compares[[n]][c(2:3)], collapse = '_vs._'),
                           '.pdf'), width=12, height = 8)
  
}

##########################################
# make heatmap  
##########################################
library(grid)
ggs = unique(ggs)
paletteLength <- 50

cpm = fpm(dds1, robust = TRUE)
cpm = cpm[match(ggs, rownames(cpm)), ]
cc = unique(design$conds)

mat = matrix(NA, nrow = nrow(cpm), ncol = length(cc))
colnames(mat) = cc
rownames(mat) = rownames(cpm)

for(n in seq_len(ncol(mat)))
{
  mat[,n] = apply(cpm[, which(design$conds == colnames(mat)[n])], 1, median)
}

mat = mat[, c(2,4,6,1,3,5)]

mat1 <- t(scale(t(mat)))

myColor <- viridis::viridis(paletteLength)
#myColor1 <- colorRampPalette(ArchRPalettes$coolwarm)(paletteLength)
#myColor2 <- colorRampPalette(c('lightgray', 'red'))(paletteLength)

myBreaks <- c(seq(min(mat1), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(mat1)/paletteLength, max(mat1), length.out=floor(paletteLength/2)))

# set annotation
anno_col <- data.frame(cbind(sortedCells =  rep(c('FoxA2+', 'FoxA2-'), times = c(3, 3)), 
                             conditions = rep(c('RA', 'BMP', 'LDN'),  time= 2)))
rownames(anno_col) <- colnames(mat1)

col<- colorRampPalette(c("blue4", "white", "darkred"))(paletteLength)
col = colorRampPalette(rev(brewer.pal(n = 7, name ="RdGy")))(paletteLength)
#col = colorRampPalette(rev(brewer.pal(n = 7, name ="PuOr")))(16)
col = palette(hcl.colors(8, "Viridis"))

col = colorRampPalette(c("navy", "white", "red3"))(paletteLength)

library(pheatmap)
pheatmap::pheatmap(
    mat1,
    border_color = 'gray60',
    color = col,
    breaks = myBreaks,
    annotation_col = anno_col,
    angle_col = 45,
    cluster_cols = FALSE,
    fontsize_row = 4,
    treeheight_row = 30,
    filename = paste0(outDir, 'heatmap_sortedCells_conditions_with.geneNames.pdf'), 
    width = 8, height = 20
)

pheatmap::pheatmap(
  mat1,
  border_color = 'gray60',
  color = col,
  breaks = myBreaks,
  annotation_col = anno_col,
  angle_col = 45,
  cluster_cols = FALSE,
  fontsize_row = 4,
  show_rownames = FALSE,
  treeheight_row = 30,
  filename = paste0(outDir, 'heatmap_sortedCells_conditions_without.geneNames.pdf'), 
  width = 5, height = 10
)

