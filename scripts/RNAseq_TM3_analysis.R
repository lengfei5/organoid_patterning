##########################################################################
##########################################################################
# Project: Organoid patterning project
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed May 12 15:13:46 2021
##########################################################################
##########################################################################
version.analysis = 'Rnaseq_old_Dresden'

resDir = paste0("../results/RNAseq.old.by.Maria.", version.analysis)
RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../data/'

design = read.delim(paste0(dataDir, 'timeSeries_SampleMetadata.txt'), header = TRUE, sep = '\t', skip = 1, comment.char = "#")
design = design[grep('*_Bigwig', design$Linking_id, invert = TRUE), ]
design = design[grep('^_', design$Linking_id, invert = TRUE), ]

cc = unique(design$condition)
length(unique(design$condition))

counts = read.delim(paste0(dataDir, 'timeSeries_featurecounts.count.gene.tsv'), header = TRUE, sep = '\t')
genes = counts$gene_id

annots = read.delim('/Volumes/groups/tanaka/People/current/jiwang/annotations/mouse/mm10/ensemble/Ensemble_biomart_mm10.txt', 
                    sep = '\t', header = TRUE)
mm = match(genes, annots$Gene.stable.ID)

genes = data.frame(c(1:nrow(counts)), genes, annots$Gene.name[mm], gene.types = annots$Gene.type[mm], stringsAsFactors = FALSE)
genes = genes[which(!is.na(genes[, 3])), ]
genes = genes[which(genes$gene.types == 'protein_coding'), ]
genes = genes[match(unique(genes[,3]), genes[,3]), ]


counts = counts[genes[,1], -1]
rownames(counts) = genes$annots.Gene.name.mm. 

design$Linking_id = sapply(design$Linking_id, function(x) {test = unlist(strsplit(as.character(x), '-')); return(test[length(test)])} )

kk = match(design$Linking_id, colnames(counts))
counts = counts[,kk]
colnames(counts) = paste0(design$condition, '_', design$Linking_id)

require(ggplot2)
require(DESeq2)
library("dplyr")
library("ggplot2")

dds <- DESeqDataSetFromMatrix(counts, DataFrame(design), design = ~ condition)

ss = rowSums(counts(dds))

hist(log2(ss), breaks = 200, main = 'log2(sum of reads for each gene)')

cutoff.peak = 2^10
cat(length(which(ss > cutoff.peak)), 'peaks selected \n')

dds <- dds[ss > cutoff.peak, ]

# normalization and dimensionality reduction
dds <- estimateSizeFactors(dds)
fpm = fpm(dds, robust = TRUE)
fpm = log2(fpm + 2^-6)
hist(fpm, breaks = 100)

# average triplicates
tt = c(-18, -10, 0, 12, 24, 36, 48, 60)

fpm.RA = matrix(NA, nrow = nrow(fpm), ncol = length(tt))
fpm.noRA = matrix(NA, nrow = nrow(fpm), ncol = length(tt))
rownames(fpm.RA) = rownames(fpm)
rownames(fpm.noRA) = rownames(fpm)

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
  
  fpm.RA[,n] = apply(fpm[,kk1], 1, median)
  fpm.noRA[,n] = apply(fpm[, kk2], 1, median)
 
}


examples = unique(c('Foxa2', # FoxA 
  'Lef1', 'Mapk1', rownames(fpm)[grep('Smad', rownames(fpm))], 
  rownames(fpm)[grep('Wnt|Dkk|Tcf', rownames(fpm))],
  rownames(fpm)[grep('Bmp', rownames(fpm))], 'Nog', 'Chrd', 'Runx1', 'Runx2', 
  rownames(fpm)[grep('Fgf', rownames(fpm))], 'Dusp1', 'Dusp10', 'Dusp27', 'Dusp4', 'Dusp5', 'Mapk10', 'Mapk4', 'Mapk8ip2', 'Spry4' 
             ))

pdfname = paste0(resDir, '/Wnt_targets_RANseq_timeSeries.pdf')
pdf(pdfname,  width = 10, height = 8)
par(cex = 1.0, las = 1, mgp = c(3,2,0), mar = c(6,6,2,0.2), tcl = -0.3)

for(g in examples)
{
 
  kk = which(rownames(fpm) == g)
  if(length(kk)){
    cat(g, '\n')
    
    plot(c(0, 1), type = 'n', xlim = c(-18, 60), ylim = range(c(fpm[kk, ], 0, 2.0, 3.5)), main = g, 
         ylab = 'log2(cpm)', xlab = 'time')
    points(tt, fpm.RA[kk, ], col = 'darkblue', type = 'l', pch = 16, lwd = 2.0)
    points(tt, fpm.RA[kk, ], col = 'darkblue', type = 'p', pch = 16)
    
    points(tt, fpm.noRA[kk, ], col = 'darkred', type = 'l', lwd = 2.0)
    points(tt, fpm.noRA[kk, ], col = 'darkred', type = 'p', pch = 1)
    
    abline(h = 2.5, col = 'darkgray', lwd = 2.0)
    legend('topright', legend = c('RA treatment', 'no treatment'), bty = 'n', 
           col = c('darkblue', 'darkred', lwd =2.0, pch = c(16, 1)), lty = 1)
    
  }else{
    cat(g, 'Not Found \n')
  }
  
}

dev.off()



