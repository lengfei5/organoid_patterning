##########################################################################
##########################################################################
# Project: Organoid patterning project
# Script purpose: analyze the time series RNA-seq data by Elena and sorted cells Foxa2+ at day 5
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed May 12 15:13:46 2021
##########################################################################
##########################################################################
rm(list = ls())
version.analysis = 'Rnaseq_old_Dresden'

resDir = paste0("../results/RNAseq.old.by.Maria.", version.analysis)
RdataDir = paste0('../results/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../data/'

require(ggplot2)
require(DESeq2)
require(dplyr)
require("pheatmap")

########################################################
########################################################
# Section : import data and metadata
# 
########################################################
########################################################
design = read.delim(paste0(dataDir, 'timeSeries_SampleMetadata.txt'), header = TRUE, sep = '\t', skip = 1, comment.char = "#")
design = design[grep('*_Bigwig', design$Linking_id, invert = TRUE), ]
design = design[grep('^_', design$Linking_id, invert = TRUE), ]

design.sorted = read.delim(paste0(dataDir, 's48h_SampleMetadata.txt'), header = TRUE, sep = '\t', skip = 1, comment.char = "#")
design.sorted = design.sorted[grep('*_Bigwig', design.sorted$Linking_id, invert = TRUE), ]
design.sorted = design.sorted[grep('^_', design.sorted$Linking_id, invert = TRUE), ]

design = rbind(design, design.sorted)

cc = unique(design$condition)
length(unique(design$condition))

counts = read.delim(paste0(dataDir, 'timeSeries_featurecounts.count.gene.tsv'), header = TRUE, sep = '\t')

counts.sorted = read.delim(paste0(dataDir, 's48h_featurecounts.count.gene.tsv'), header = TRUE, sep = '\t')
ggs = counts.sorted$gene_id
ggs = sapply(ggs, function(x) unlist(strsplit(as.character(x), '[|]'))[1])
counts.sorted$gene_id = ggs

mm = match(counts$gene_id, counts.sorted$gene_id)
counts = cbind(counts, counts.sorted[mm, -1])

rm(design.sorted)
rm(counts.sorted)

save(counts, design, file = paste0(RdataDir, '/RNAseq_timeSeries_sortedDay5_count_design.Rdata'))

##########################################
# annnotation: convert ens ID to gene symbol 
##########################################
load(file = paste0(RdataDir, '/RNAseq_timeSeries_sortedDay5_count_design.Rdata'))
genes = counts$gene_id

annots = read.delim('/Volumes/groups/tanaka/People/current/jiwang/Genomes/mouse/mm10_ens/ens_BioMart_GRCm38.p6.txt', sep = '\t', 
                    header = TRUE)
annots = annots[which(annots$Gene.type == 'protein_coding'), ]

mm = match(genes, annots$Gene.stable.ID)

# extract gene lengths by calculating the median of transcript length
Find.gene.length = FALSE
if(Find.gene.length){
  genes = data.frame(c(1:nrow(counts)), genes, annots$Gene.name[mm], gene.types = annots$Gene.type[mm], stringsAsFactors = FALSE)
  genes = genes[which(!is.na(genes[, 3])), ]
  
  genes = genes[match(unique(genes[,3]), genes[,3]), ]
  colnames(genes) = c('index.counts', 'ensID', 'gene', 'genetype')
  
  genes$length = NA
  for(n in 1:nrow(genes))
  {
    if(n%%100 ==0) cat(n, '\n')
    genes$length[n] = median(annots$Transcript.length..including.UTRs.and.CDS.[which(annots$Gene.name == genes$gene[n])])
  }
  
  saveRDS(genes, file = paste0(RdataDir, 'RNAseq_timeSeries_Foxa2.pos.day5_gene_names_length_.rds'))
  
}else{
  genes = readRDS(file = paste0(RdataDir, '/gene_names_length.rds'))
}

counts = counts[genes[,1], -1]
rownames(counts) = genes$gene

design$Linking_id = sapply(design$Linking_id, function(x) {test = unlist(strsplit(as.character(x), '-')); return(test[length(test)])} )

colnames(counts)[46:ncol(counts)] = sapply(colnames(counts)[46:ncol(counts)], function(x) unlist(strsplit(as.character(x), '[.]'))[2])
kk = match(design$Linking_id, colnames(counts))
counts = counts[,kk]
colnames(counts) = paste0(design$condition, '_', design$Linking_id)

save(counts, design, genes, file = paste0(RdataDir, '/RNAseq_timeSeries_sortedDay5_count_design_geneSymbol.Rdata'))

########################################################
########################################################
# Section :  DESeq2 normalization
# 
########################################################
########################################################
load(file = paste0(RdataDir, '/RNAseq_timeSeries_sortedDay5_count_design_geneSymbol.Rdata'))

dds <- DESeqDataSetFromMatrix(counts, DataFrame(design), design = ~ condition)

ss = rowSums(counts(dds))

hist(log2(ss), breaks = 200, main = 'log2(sum of reads for each gene)')

cutoff.peak = 100
cat(length(which(ss > cutoff.peak)), 'peaks selected \n')
gg.tokeep = 'Nog|Acvrl1|Acvr1|Notch|Dll|Jagn|Gdf|Bmp|Fgf|Wnt'

sels = unique(c(which(ss > cutoff.peak), grep(gg.tokeep, rownames(dds))))

dds <- dds[sels, ]

# normalization and dimensionality reduction
dds <- estimateSizeFactors(dds)
fpm = fpm(dds, robust = TRUE)
ss = colSums(counts(dds))
plot(sizeFactors(dds), ss/10^6)

#dds0 = dds
#fpm0 = fpm
#design0 = design
#save(design0, dds0, fpm0, file = paste0(RdataDir, '/RNAseqOld_design_dds_fpm.Rdata'))
save(design, dds, file = paste0(RdataDir, '/RNAseq_timeSeries_sortedDay5_count_design_geneSymbol_dds.Rdata'))

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

#kk = grep('RA', design$condition)
kk = c(1:nrow(design))
pca=plotPCA(vsd[,kk], intgroup = c('condition'), returnData = TRUE, ntop = 500)
pca2save = as.data.frame(pca)
ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition))  + 
  geom_point(size=3) + 
  geom_text(hjust = 0.3, nudge_y = 0.4, size=4)

plot(ggp) + ggsave(paste0(resDir, "/PCAplot_RA.noRA.timeSeries_Foxa2pos.pdf"), width=18, height = 12)


kk = grep('RA', design$condition)
#kk = c(1:nrow(design))
pca=plotPCA(vsd[,kk], intgroup = c('condition'), returnData = TRUE, ntop = 500)
pca2save = as.data.frame(pca)
ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition))  + 
  geom_point(size=3) + 
  geom_text(hjust = 0.3, nudge_y = 0.4, size=4)

plot(ggp) + ggsave(paste0(resDir, "/PCAplot_RA.timeSeries_Foxa2pos.pdf"), width=18, height = 12)

########################################################
########################################################
# Section : genome-wide anaysis to identify pattern-related genes from time-series RNA-seq data 
# 
########################################################
########################################################
load(file = paste0(RdataDir, '/RNAseq_timeSeries_sortedDay5_count_design_geneSymbol_dds.Rdata'))

fpm = fpm(dds, robust = TRUE)

## RA treatment only
kk = setdiff(grep('RA', design$condition), grep('RA_AF|RA_GFPp', design$condition))
dds1 = dds[,kk]
dds1$condition = droplevels(dds1$condition)
res1 = DESeq(dds1, test="LRT", reduced=~1)

res1 <- as.data.frame(results(res1))

fpm1 = fpm[, kk]
xx = design[kk, ]
cpm = matrix(NA, nrow = nrow(fpm1), ncol = length(unique(xx$condition)))
rownames(cpm) = rownames(fpm1)
cc = unique(xx$condition)
cc = cc[c(1,5,6,2,7,8,3,4)]
colnames(cpm) = cc
for(n in 1:ncol(cpm))
{
  # n = 1
  jj = which(xx$condition == colnames(cpm)[n])
  cpm[,n] = apply(fpm1[,jj], 1, median)
}

## RA treatment and control
kk = grep('s48h|before_RA', design$condition, invert = TRUE)
dds2 = dds[,kk]
xx = design[kk,]
xx$condition = as.character(xx$condition)
xx$time = xx$condition

xx$condition[grep('RA', xx$condition, invert = TRUE)] = 'control'
xx$condition[grep('RA', xx$condition)] = 'RA'
xx$time = gsub('_RA', '', xx$time)

ddsTC <- DESeqDataSetFromMatrix(counts(dds2), data.frame(xx), ~ condition + time + condition:time)
sizeFactors(ddsTC) = sizeFactors(dds2)

ddsTC <- DESeq(ddsTC, test="LRT", reduced = ~ condition + time)
resTC <- results(ddsTC)

save(cpm, res1, resTC, file = paste0(RdataDir, '/RNAseq_timeSeries_fpmMean_DEtestRes.timepoints.and.vs.control.Rdata'))

## combine genes with non-static padj < 0.01 and behaving different from control condition padj <0.01
jj = which(res1$padj< 0.01 & resTC$padj< 0.01)

yy = log2(cpm[jj,] + 2^-6)
diff = apply(yy, 1, max) - apply(yy, 1, min)
yy = yy[which(diff>=1), ]

df <- data.frame(cond = colnames(yy), time = c(rep('day2', 3), rep('day3', 2), rep('day4', 2), rep('day5', 1)))
rownames(df) = colnames(yy)

pheatmap(yy, cluster_rows=TRUE, show_rownames=FALSE, scale = 'row', show_colnames = TRUE, 
         cluster_cols=FALSE, na_col = 'gray',  annotation_col = df, gaps_col = c(3, 5, 7),
         filename = paste0(resDir, '/timeSeries_variableGenes.pdf'), width = 10, height = 8)

##########################################
# pathways analysis 
##########################################
library(org.Mm.eg.db)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(stringr)

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

gg.expressed = rownames(yy)
  
bgs = readRDS(file = paste0(RdataDir, 'RNAseq_timeSeries_Foxa2.pos.day5_gene_names_length_.rds'))
bgs = bgs$gene
gene.df <- bitr(gg.expressed, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Mm.eg.db)
head(gene.df)

bgs.df <- bitr(bgs, fromType = "SYMBOL",
               toType = c("ENSEMBL", "ENTREZID"),
               OrgDb = org.Mm.eg.db)
head(bgs.df)


ego <-  enrichGO(gene         = gene.df$ENSEMBL,
                 universe     = bgs.df$ENSEMBL,
                 #OrgDb         = org.Hs.eg.db,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'ENSEMBL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05, readable = TRUE)

#head(ego)
barplot(ego, showCategory=50) + ggtitle("Go term enrichment for time-series RNA-seq data") +  
  ggsave(paste0(resDir, "/RNAseq_timeSeries_siganificantGenes.dynamicGenes.logFC.1_enrichGO_top50.pdf"), width=12, height = 12)
  
ii = grep('pathway', ego[, 2])
xx = ego[ii, ]
xx = xx[order(-xx$p.adjust), ]
xx$Description = gsub('signaling', 'S', xx$Description)
xx$Description = gsub(' pathway', 'P', xx$Description)
par(cex = 1.0, las = 1, mgp = c(3,1,0), mar = c(3, 30,2,0.2), tcl = -0.3)
barplot(xx$Count,  names.arg = xx$Description, horiz = TRUE, beside = TRUE, las=2)

write.table(data.frame(ego), 
            file = paste0(resDir, '/enrichGO_RNAseq_timeSeries_siganificantGenes.dynamicGenes.logFC.1_enrichGO_top50.txt'),
            sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)

write.table(xx, 
            file = paste0(resDir, '/enrichGO_RNAseq_timeSeries_siganificantGenes.dynamicGenes.logFC.1_singalingPathways.txt'),
            sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
# kegg enrich did not work quite well
# kegg.up <- enrichKEGG(gene         =  gene.df$ENTREZID,
#                  organism     = 'mmu',
#                  universe     = bgs.df$ENTREZID,
#                  keyType = 'kegg',
#                  pvalueCutoff = 0.05)
# head(kegg.up)
# barplot(kegg.up, showCategory=50) + ggtitle("kegg for upregulated genes")

Test.signaling.pathways.for.each.timepoint = FALSE
if(Test.signaling.pathways.for.each.timepoint){
  load(file = paste0(RdataDir, '/RNAseq_timeSeries_sortedDay5_count_design_geneSymbol_dds.Rdata'))
  
  kk = grep('s48h', design$condition, invert = TRUE)
  dds = dds[,kk]
  dds$condition = droplevels(dds$condition)
  
  dds <- DESeq(dds)
  resultsNames(dds)  
  
  cc = unique(dds$condition)  
  cc = cc[grep('RA', cc)]
  cc = cc[c(1, 5, 6, 2, 7, 8,3, 4)]
  
  for(n in 2:length(cc))
  {
    # n = 2
    cat(as.character(cc[n]), '\n')
    res = results(dds, contrast=c("condition", as.character(cc[n]), as.character(cc[n-1])))
    ggs = rownames(res)[which(res$padj < 0.05 & res$log2FoldChange > 0.5 )]
    gene.df <- bitr(ggs, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Mm.eg.db)
    ego <-  enrichGO(gene         = gene.df$ENSEMBL,
                     universe     = bgs.df$ENSEMBL,
                     #OrgDb         = org.Hs.eg.db,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = 'ENSEMBL',
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05, readable = TRUE)
    
    #barplot(ego, showCategory=50) 
    ii = grep('pathway', ego[, 2])
    xx = ego[ii, ]
    
    write.table(xx, 
                file = paste0(resDir, '/enrichGO_RNAseq_timepoint_', cc[n], '.DEgenes.logFC.0.5_singalingPathways.txt'),
                sep = '\t', col.names = TRUE, row.names = TRUE, quote = FALSE)
    
  }
  
  files = list.files(path = resDir, pattern = 'RA.DEgenes.logFC.0.5_singalingPathways.txt', full.names = TRUE)
  files = files[c(7, 1:6)]
  
  ids = c()
  for(n in 1:length(files))
  {
    # n = 1
    ff = read.table(files[n], sep = '\t', header = TRUE, row.names = 1, as.is = c(1, 2,8))
    ids = rbind(ids, cbind(ff$ID, as.character(ff$Description), as.character(ff$geneID) ))
  }
  
  ids = data.frame(ids)
  colnames(ids) = c('ID', 'description', 'geneID')
  
  ##########################################
  # keep the relevant significant genes for signaling pathways 
  ##########################################
  id.uniq = unique(ids$ID)
  xx = ids[match(id.uniq, ids$ID), ]
  xx$geneID = as.character(xx$geneID)
  for(n in 1:nrow(xx))
  {
    # n = 1
    jj = which(ids$ID == xx$ID[n])
    if(length(jj)>1) {
      gg = c()
      for(j in jj) gg = c(gg, unlist(strsplit(as.character(ids$geneID[j]), '/')))
      gg = unique(gg)
      #gg = sapply(ids$geneID[jj], function)
      xx$geneID[n] = paste0(gg, collapse = ',')
      
    }
  }
  
  ids = xx
  
  res = matrix(NA, nrow = nrow(ids), ncol = length(files))
  rownames(res) = ids$ID
  tt = c()
  for(n in 1:length(files))
  {
    # n = 1
    ff = read.table(files[n], sep = '\t', header = TRUE, row.names = 1, as.is = c(1, 2,8))
    kk = match(ff$ID, rownames(res))
    res[kk, n] = ff$pvalue
    tt = c(tt, gsub('_RA.DEgenes.logFC.0.5_singalingPathways.txt', '', gsub('enrichGO_RNAseq_timepoint_', '', basename(files[n]))))
  }
  res = data.frame(res)
  colnames(res) = tt
  mm = match(rownames(res), ids$ID)
  rownames(res) = ids$description[mm]
  
  yy = -log10(res)
  pheatmap(yy, cluster_rows=FALSE, show_rownames=TRUE, scale = 'none', show_colnames = TRUE, 
           cluster_cols=FALSE, na_col = 'gray', gaps_col = c(2, 4, 6), fontsize_row = 12,
           filename = paste0(resDir, '/siganling_pathways_activationTiming_stringent_logFC0.5_FDR0.01.pdf'), width = 16, height = 12)
  
}

########################################################
########################################################
# Section : Normalized data with gene length and plot examples
# 
########################################################
########################################################
load(file = paste0(RdataDir, '/RNAseq_timeSeries_sortedDay5_count_design_geneSymbol_dds.Rdata'))
genes = readRDS(file = paste0(RdataDir, '/RNAseq_timeSeries_Foxa2.pos.day5_gene_names_length.rds'))
load(file = paste0(RdataDir, '/RNAseq_timeSeries_fpmMean_DEtestRes.timepoints.and.vs.control.Rdata'))
#fpm = fpm(dds, robust = TRUE)
#sps = read.table(file = paste0(resDir, '/enrichGO_RNAseq_timeSeries_siganificantGenes.dynamicGenes.logFC.1_singalingPathways.txt'),
#sep = '\t', header = TRUE, row.names = 1)
ggs = readRDS(file = paste0(RdataDir, '/TM3_examplesGenes_withGOterm.rds'))

##########################################
# normalize with gene length 
##########################################
cpm = fpm(dds, robust = TRUE)
rpkm = cpm
ll = genes$length[match(rownames(cpm), genes$gene)]
for(n in 1:ncol(rpkm))
{
  rpkm[,n] = cpm[,n] / ll *10^3
}


rpkm = log2(rpkm + 2^-6)

hist(rpkm, breaks = 200)
abline(v = c(1, 2), col = 'red', lwd = 2.0)

# average triplicates0
tt = c(-18, -10, 0, 12, 24, 36, 48, 60)

rpkm.RA = matrix(NA, nrow = nrow(rpkm), ncol = length(tt))
rpkm.noRA = matrix(NA, nrow = nrow(rpkm), ncol = length(tt))
rownames(rpkm.RA) = rownames(rpkm)
rownames(rpkm.noRA) = rownames(rpkm)

for(n in 1:length(tt))
{
  if(tt[n] == -18){
    kk1 = which(design$condition == 'before_RA')
    kk2 = which(design$condition == 'before_RA')
  }else{
    kk2 = grep(paste0(as.character(abs(tt[n]), 'h')), design$condition)
    kk1 = kk2[grep('RA', design$condition[kk2])]
    kk2 = kk2[grep('RA', design$condition[kk2], invert = TRUE)]
  }
  
  rpkm.RA[,n] = apply(rpkm[,kk1], 1, median)
  rpkm.noRA[,n] = apply(rpkm[, kk2], 1, median)
}


##########################################
# select only signaficant changes genes and normalized with transcript length
# make a summary plots with potential players from Wnt, FGF, BMP
##########################################
Select.signficant.Genes = FALSE
if(Select.signficant.Genes){
  # rpkm = cpm
  # ll = genes$length[match(rownames(cpm), genes$gene)]
  # for(n in 1:ncol(rpkm))
  # {
  #   rpkm[,n] = cpm[,n] / ll *10^3
  # }
  # 
  # rpkm = log2(rpkm + 2^-6)
  # hist(rpkm, breaks = 200)
  # abline(v = c(1, 2), col = 'red', lwd = 2.0)
  # 
  
  ### we will select genes: showing changes across time points and different from noRA treatment
  vars = apply(rpkm.RA, 1, var)
  means = apply(rpkm.RA, 1, mean)
  df = rpkm.RA - rpkm.noRA
  df = apply(df[, -1], 1, mean)
  
  plot(df, vars)
  
  xx = data.frame(gene = names(df), diff = df, vars = vars, stringsAsFactors = FALSE)
  
  mm = match(xx$gene, examples)
  xx = xx[!is.na(mm), ]
  
  pdfname = paste0(resDir, "/genes_signalingPathways_variances_RA.noRA.diff.pdf")
  pdf(pdfname, width = 16, height = 12)
  
  p0 = ggplot(xx, aes(x = diff, y = vars, label = gene)) +
    geom_point(size = 0.5) + 
    geom_text(data=subset(xx, diff > 2. | diff < -2. | vars > 1.5), size = 4, nudge_y = 0.1) + 
    theme(legend.position = "right") 
  plot(p0) 
    #geom_text(aes(label = gg.examepls, x = 12, colour = names(lab), y = c(lab), hjust = -.02))
  dev.off()
  
  jj = which(res1$padj< 0.01 & resTC$padj< 0.01)
  rpkm = rpkm[jj, ]
  
  Make.summary.plot = FALSE
  if(Make.summary.plot){
    require(ggplot2)
    require(grid)
    require(gridExtra)
    library(tidyr)
    library(dplyr)
    library(ggrepel)
    library(RColorBrewer)
    
    tt = c(-18, -10, 0, 12, 24, 36, 48, 60)
    
    #gg.examepls = c('Foxa2', 'Lef1','Wnt3', 'Sfrp5',  'Fgf2', 'Fgf8', 'Bmp7', 'Bmp4', 'Bmp6', 'Nog')
    mains = 'positive_controls'
    gg.examepls = c('Foxa2', 'Shh', 'Olig2')
    #gg.examepls = c('Lef1', 'Wnt1', 'Wnt3', 'Wnt3a', 'Wnt7a', 'Wnt7b', 'Wnt8a') 
    #gg.examepls = c('Lypd6', 'Dkk3', 'Sfrp5', 'Dkk1', 'Sost')
    #gg.examepls = c('Spry4', 'Spry2', 'Etv4', 'Etv5')
    #gg.examepls = c('Fgf10', 'Fgf17','Fgf8', 'Fgf5', 'Fgf2','Fgf21', 'Fgf11','Fgf1', 'Fgf4')
    gg.examepls = c('Id1', 'Id3', 'Smad6', 'Nog', 'Fst', 'Bambi')
    #gg.examepls = c('Bmp7', 'Bmp4', 'Bmp1', 'Bmp6', 'Bmpr2', 'Bmpr1b')
    
    pdfname = paste0(resDir, "/genes_signalingPathways_summaryPlots_Foxa2.pdf")
    pdf(pdfname, width = 12, height = 8)
    par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
    
    cols = colorRampPalette(rev(brewer.pal(9, "RdBu")))(length(gg.examepls)+1)
    xx = rpkm.RA[match(gg.examepls, rownames(rpkm.RA)), ]
    #xx = data.frame(xx, stringsAsFactors = FALSE)
    #df = as_tibble(xx) %>% gather(time, expression, 3:10)
    
    plot(c(0, 1), type = 'n', xlim = c(-22, 62), ylim = range(c(xx, 0, 1), na.rm = TRUE), 
         ylab = 'log2(RPKM)', xlab = 'time', main = mains)
    for(kk in 1:nrow(xx)){
      points(tt, xx[kk, ], type = 'l', pch = 16, lwd = 2.0, col = cols[which(gg.examepls == gg.examepls[kk])])
      points(tt, xx[kk, ], type = 'p', pch = 16, cex = 1.2,col = cols[which(gg.examepls == gg.examepls[kk])])
      text(tt[1], xx[kk, 1], labels = rownames(xx)[kk], pos = 2, offset = 0.5)
    }
    abline(h = c(0, 1), col = 'darkgray', lwd = 1.5, lty = 3)
    abline(v = c(6, 32, 54), col = 'cornflowerblue', lwd = 2.0, lty=3)
    
    dev.off()
    
    # ggplot(df, aes(x = time, y = expression, group = gene, color = pathways)) +
    #   geom_line(size=1.2, linetype = 'solid') +
    #   geom_point(size = 2.5) + 
    #   theme(legend.position = "right") +
    #   geom_text(aes(label = gg.examepls, x = 12, colour = names(lab), y = c(lab), hjust = -.02))
    # 
    #   geom_bar(stat = "identity") +
    #   theme(legend.position = "none")  + 
    #   ggtitle('nb of cysts ') + 
    #   theme(axis.text.x = element_text(angle = 90))
    
  }
  
  MakeHeatmap = FALSE
  if(MakeHeatmap){
    diff = apply(rpkm, 1, max) - apply(rpkm, 1, min)
    hist(diff, breaks = 100)
    abline(v = c(1, 0.5), col = 'red', lwd = 2.0)
    length(which(diff>1))
    length(which(diff>0.5))
    
    rpkm = rpkm[which(diff>=1), ]
    
    kk = match(rownames(rpkm), ggs$gene)
    rpkm = rpkm[which(!is.na(kk)), ]
    
    rpkm = rpkm[c(which(rownames(rpkm) == 'Foxa2'), 
                  which(!is.na(match(rownames(rpkm), ggs$gene[which(ggs$pathway == 'WNT')]))),
                  which(!is.na(match(rownames(rpkm), ggs$gene[which(ggs$pathway == 'FGF')]))),
                  which(!is.na(match(rownames(rpkm), ggs$gene[which(ggs$pathway == 'BMP')]))),
                  which(!is.na(match(rownames(rpkm), ggs$gene[which(ggs$pathway == 'NOTCH')]))),
                  which(!is.na(match(rownames(rpkm), ggs$gene[which(ggs$pathway == 'SHH')])))
    ), ]
    
    yy = rpkm
    max = apply(yy, 1, max)
    yy = yy[which(max>1), ]
    #yy[which(yy<0)] = NA
    pheatmap(yy, cluster_rows=TRUE, show_rownames=TRUE, scale = 'row', show_colnames = TRUE, 
             cluster_cols=FALSE, na_col = 'gray', gaps_col = c(3, 5, 7), fontsize_row = 8,
             filename = paste0(resDir, '/siganling_pathways_geneExpression_timeSeries.pdf'), 
             width = 12, height = 28)
    
  }
  
 
}

##########################################
# make example plots of WNT, BMP, FGF ligand, receptor, effectors and targets
##########################################
Make.example.plots = FALSE
if(Make.example.plots){
  
  sorted = cbind(apply(rpkm[ , grep('s48h_RA_AF', colnames(rpkm))], 1, median), 
                 apply(rpkm[ , grep('s48h_RA_GFPp', colnames(rpkm))], 1, median),
                 apply(rpkm[ , grep('s48h_SAG_GFPp', colnames(rpkm))], 1, median))
  
  
  examples = unique(c('Foxa2', # FoxA 
                      'Lef1', 'Mapk1', rownames(rpkm)[grep('Smad', rownames(rpkm))], 
                      rownames(rpkm)[grep('Wnt|Dkk|Tcf|Axin|Ctnnb|Sfr|Rspo|Igfbp|Wif|Ndp', rownames(rpkm))],
                      rownames(rpkm)[grep('Bmp|Gdf', rownames(rpkm))], 'Nog', 'Chrd', 'Runx1', 'Runx2',  'Smad6', 'Id1', 'Id3',
                      rownames(rpkm)[grep('Acvr', rownames(rpkm))],
                      rownames(rpkm)[grep('Fgf|Etv|Il17rd', rownames(rpkm))], 
                      'Dusp1', 'Dusp10', 'Dusp27', 'Dusp4', 'Dusp5', 'Mapk10', 'Mapk4', 'Mapk8ip2', 'Spry4', 'Rbpj', 'Hes1', 'Hes5',
                      'Hes7', 'Hey1', 'Hey2',
                      rownames(rpkm)[grep('Notch|Jag|Dll|Dlk', rownames(rpkm))], 
                      rownames(rpkm)[grep('Tgf', rownames(rpkm))]
  ))
  
  examples = unique(c(examples, ggs$gene))
  
  cat('nb of genes to check :', length(examples), '\n')
  
  pdfname = paste0(resDir, '/RANseq_timeSeries_sortedFoxA2positive_genes_pathways_v7.pdf')
  pdf(pdfname,  width = 10, height = 6)
  par(cex = 1.0, las = 1, mgp = c(3,2,0), mar = c(6,6,2,0.2), tcl = -0.3)
  
  hist(rpkm.RA, breaks = 100, xlab = 'log2(RPKM)')
  abline(v = c(0, 1), col = 'red', lwd = 2.0)
  
  hist(rpkm.noRA, breaks = 100, xlab = 'log2(RPKM)')
  abline(v = c(0, 1), col = 'red', lwd = 2.0)
  
  hist(sorted, breaks = 100, xlab = 'log2(RPKM)')
  abline(v = c(0, 1), col = 'red', lwd = 2.0)
  
  
  for(g in examples)
  {
    # g = 'Foxa2'
    kk = which(rownames(rpkm) == g)
    if(length(kk)){
      cat(g, '\n')
      
      plot(c(0, 1), type = 'n', xlim = c(-18, 60), ylim = range(c(rpkm[kk, ], 0, 1)), main = g, 
           ylab = 'log2(RPKM)', xlab = 'time')
      points(tt, rpkm.RA[kk, ], col = 'darkblue', type = 'l', pch = 16, lwd = 2.0)
      points(tt, rpkm.RA[kk, ], col = 'darkblue', type = 'p', pch = 16)
      
      points(tt, rpkm.noRA[kk, ], col = 'darkred', type = 'l', lwd = 1.0)
      points(tt, rpkm.noRA[kk, ], col = 'darkred', type = 'p', pch = 1)
      
      points(48, sorted[kk, 1], col = 'darkblue', type = 'p', cex = 2.0, pch = 21, bg = 'magenta')
      points(48, sorted[kk, 2], col = 'darkblue', type = 'p', cex = 2.0, pch = 21, bg = 'darkgreen')
      #points(48, sorted[kk, 3], col = 'darkorange', type = 'p', cex = 2.0, pch = 18)
      
      abline(h = c(0, 1), col = 'darkgray', lwd = 2.0)
      abline(v = c(6, 32, 54), col = 'cornflowerblue', lwd = 2.0, lty=3)
      legend('topleft', legend = c('RA', 'no.RA', 's48h.RA.AF', 's48h.RA.GFPp'), bty = 'n', 
             col = c('darkblue', 'darkred', 'magenta', 'darkgreen'), lwd =2.0, pch = c(16, 1, 16, 16), lty = c(1, 1, 0, 0))
      
    }else{
      cat(g, 'Not Found \n')
    }
  }
  
  dev.off()
  
}




