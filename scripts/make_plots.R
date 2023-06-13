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

Counts.to.Use = 'UMI'

load(file = paste0(RdataDir, 'TM3_dds_normalized.Rdata'))
#load(file = paste0(resDir, '/TM3_Foxa2.positive_vs_neg.day3_comparisonStats_cpm.Rdata'))


load(file = paste0(RdataDir, '/design.detailed_RawUMI_', Counts.to.Use, version.analysis, '.Rdata'))

load(file = paste0(RdataDir, '/TM3_pooled.pos.neg_ddx_cc.pools_', Counts.to.Use, version.analysis, '.Rdata'))

cc = unique(design.matrix$condition)
cc = cc[which(cc != "N2B27" & cc != 'RA')]


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






