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
require(ggplot2)
require(DESeq2)
require(dplyr)

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

dds0 = dds
fpm0 = fpm
design0 = design

save(design0, dds0, fpm0, file = paste0(RdataDir, '/RNAseqOld_design_dds_fpm.Rdata'))

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)


#kk = grep('RA', design$condition)
kk = c(1:nrow(design))
pca=plotPCA(vsd[,kk], intgroup = c('condition'), returnData = TRUE, ntop = 500)
pca2save = as.data.frame(pca)
ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition))  + 
  geom_point(size=3) + 
  geom_text(hjust = 0.3, nudge_y = 0.4, size=4)

plot(ggp) + ggsave(paste0(resDir, "/PCAplot_RA.noRA.timeSeries_Foxa2pos.pdf"), width=18, height = 12)

##########################################
# normalized by gene length
##########################################
ll = genes$length[match(rownames(fpm), genes$gene)]

rpkm = fpm

for(n in 1:ncol(rpkm))
{
  rpkm[,n] = fpm[,n] / ll *10^3
}

rpkm = log2(rpkm + 2^-6)

hist(rpkm, breaks = 100)
abline(v = 1, col = 'red', lwd = 2.0)

# average triplicates
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


sorted = cbind(apply(rpkm[ , grep('s48h_RA_AF', colnames(rpkm))], 1, median), 
               apply(rpkm[ , grep('s48h_RA_GFPp', colnames(rpkm))], 1, median),
               apply(rpkm[ , grep('s48h_SAG_GFPp', colnames(rpkm))], 1, median))


##########################################
# make plots of WNT, BMP, FGF ligand, receptor, effectors and targets
##########################################
examples = unique(c('Foxa2', # FoxA 
                    'Lef1', 'Mapk1', rownames(rpkm)[grep('Smad', rownames(rpkm))], 
                    rownames(rpkm)[grep('Wnt|Dkk|Tcf|Axin|Ctnnb', rownames(rpkm))],
                    rownames(rpkm)[grep('Bmp|Gdf', rownames(rpkm))], 'Nog', 'Chrd', 'Runx1', 'Runx2',  'Smad6', 'Id1', 'Id3',
                    rownames(rpkm)[grep('Acvr', rownames(rpkm))],
                    rownames(rpkm)[grep('Fgf', rownames(rpkm))], 
                    'Dusp1', 'Dusp10', 'Dusp27', 'Dusp4', 'Dusp5', 'Mapk10', 'Mapk4', 'Mapk8ip2', 'Spry4', 'Rbpj', 'Hes1', 'Hes5',
                    'Hes7', 'Hey1', 'Hey2',
                    rownames(rpkm)[grep('Notch|Jag|Dll|Dlk', rownames(rpkm))], 
                    rownames(rpkm)[grep('Tgf', rownames(rpkm))]
))

pdfname = paste0(resDir, '/RANseq_timeSeries_sortedFoxA2positive_genes_pathways_v5.pdf')
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
