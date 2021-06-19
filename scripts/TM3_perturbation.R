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


# setup for data import and sequencing QCs
version.analysis = '_R11601_TM3_20210619'

resDir = paste0("../results/RNAseq_perturbation", version.analysis)
RdataDir = paste0(resDir, '/Rdata')
if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/patterning_organoid/R11601/'

##########################################
# prepare design matrix, count table and QC table
##########################################
design = read.csv(paste0(dataDir, 'sampleInfo.csv'), header = TRUE)
design$condition = as.character(design$condition)
design$cells = design$condition

design$cells[grep('[+]', design$cells)] = 'Foxa2+'
design$cells[grep('[-]', design$cells)] = 'Foxa2-'
design$condition = sapply(design$condition, function(x) unlist(strsplit(as.character(x), '[_]'))[1])


stats = read.delim(paste0(dataDir, 'nf_out_RNAseq/MultiQC/multiqc_data/multiqc_general_stats.txt'), sep = '\t', 
                   header = TRUE)
alignment = read.delim(paste0(dataDir, 'nf_out_RNAseq/MultiQC/multiqc_data/multiqc_hisat2.txt')
                       , sep = '\t')

stats = stats[, c(1, 2, 3, 4, 6, 8, 9, 10)]
colnames(stats) = c('sample', 'pct.duplication', 'pct.GC', 'avg.seq.length', 'total.reads', 
                    'pct.assign', 'assigned.reads', 'alignment.rate')

stats = data.frame(stats, alignment[match(stats$sample, alignment$Sample), c(2, 4, 5)], stringsAsFactors = FALSE)
colnames(stats)[c(9:11)] = c('trimmed.reads', 'unique.aligned', 'multimapper')
stats = stats[, c(1:5, 9, 8, 10, 11, 7)]

#design = design[order(design$fileName), ]
ii = c()
jj = c()
for(n in 1:nrow(design))
{
  # n = 1;
  cat(n, '\n')
  kk = grep(design$sampleID[n], stats$sample)
  ii = c(ii, rep(n, length(kk)))
  jj = c(jj, kk)
  #kk = c(kk, grep(design$sampleID[n], stats$sample))
  #kk = c(kk, which(design$sampleID == ))
}

xx = data.frame(design[ii, ], stats[jj, ], stringsAsFactors = FALSE)
xx = xx[order(xx$fileName), ]

write.csv(xx, file = paste0(resDir, '/QCs_stats.csv'), row.names = FALSE)

design = data.frame(design, stats[kk, ], stringsAsFactors = FALSE)


########################################################
# saturation curve from rseqc
# the r code from rseqc output
########################################################
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


##################################################
## Import design matrix and prepare count table
##################################################
Processing.design.matrix = TRUE
if(Processing.design.matrix){
  library(openxlsx)
  design = read.xlsx(design.file, sheet = 1)
}

colnames(design)[1] = 'SampleID'
design$conds = paste0(design$conditions, '_', design$SampleID,  '.batch', design$batch)

# prepare the data table for different batches
batch1 = read.delim('../Data/R10724_rnaseq/nf_out_RNAseq/featureCounts/merged_gene_counts.txt', sep = '\t',  header = TRUE)
batch2 = read.delim('../Data/Rxxxx_rnaseq_old/nf_out_RNAseq/featureCounts/merged_gene_counts.txt', sep = '\t',  header = TRUE)

batch3 = batch2[, grep('HLVGMDRXX', colnames(batch2), invert = TRUE)]
batch2 = batch2[, c(1, grep('HLVGMDRXX', colnames(batch2)))]

source(RNA.functions)
xx1 = process.countTable(all=batch1, design = design[which(design$batch == '3'), c(1:2)], ensToGeneSymbol = FALSE)
colnames(xx1)[-1] = paste0(colnames(xx1)[-1], '.batch3')

xx2 = process.countTable(all=batch2, design = design[which(design$batch == '2'), c(1:2)], merge.technicalRep.sameID = TRUE,
                         ensToGeneSymbol = FALSE)
colnames(xx2)[-1] = paste0(colnames(xx2)[-1], '.batch2')

xx3 = process.countTable(all=batch3, design = design[which(design$batch == '1'), c(1:2)], merge.technicalRep.sameID = TRUE,
                         ensToGeneSymbol = FALSE)

colnames(xx3)[-1] = paste0(colnames(xx3)[-1], '.batch1')

all <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), list(xx1, xx2, xx3))

# make sure that the sample order is the same in design and all matrix
mm = match(design$conds, colnames(all))
all = all[, c(1, mm)]

write.csv(design, file = paste0(resDir, '/designSampleInfos_QCstat.csv'), row.names = FALSE)

save(design, all, file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))

##########################################
# add stats into design 
##########################################
load(file=paste0(RdataDir, 'Design_Raw_readCounts_', version.analysis, '.Rdata'))
stat0 = read.csv(file = '../results/R10724_rnaseq_202102019/QCs_stats.csv')
stat1 = read.csv(file = '../results/Rxxxx.rnaseq.old_202102020/QCs_stats.csv')

xx = data.frame(design, matrix(NA, nrow = nrow(design), ncol = 9), stringsAsFactors = FALSE)
colnames(xx)[7:15] = colnames(stat0)[-c(1:3)]

mm = match(stat0$sampleID, design$SampleID)

xx[mm, c(7:15)] = stat0[, -c(1:3)]

for(n in 1:nrow(xx))
{
  if(is.na(xx$pct.duplication[n])){
    kk = which(stat1$sampleID == xx$SampleID[n])
    if(xx$batch[n] == 1) {
      kk = kk[grep('HLVGMDRXX_', stat1$sample[kk], invert = TRUE)]
      if(length(kk) != 1) cat(n, ' -- Error \n')
      xx[n, c(7:15)] = stat1[kk, c(4:12)]
    }
    if(xx$batch[n] == 2) {
      kk = kk[grep('HLVGMDRXX_', stat1$sample[kk])]
      if(length(kk) != 2) cat(n, ' -- Error \n')
      xx[n, c(7:9, 12)] = apply(stat1[kk, c(4:6, 9)], 2, mean)
      xx[n, c(10:11, 13:15)] = apply(stat1[kk, c(7:8, 10:12)], 2, sum)
    }
    
  }
}

xx = data.frame(xx, stringsAsFactors = FALSE)

xx$pct.assigned.to.features = xx$assigned.reads/xx$total.reads *100
xx$alignment.uniq.rate = xx$unique.aligned/xx$trimmed.reads *100
xx = xx[, c(1:11, 13:15, 12, 17, 16)]

design = xx


save(design, all, file=paste0(RdataDir, 'Design_stats_readCounts_', version.analysis, '.Rdata'))

########################################################
########################################################
# Section II : analyze the RNA-seq data 
# 
########################################################
########################################################
rm(list = ls())

RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

version.Data = 'rnaseq_RNAseqSamples_all';
version.analysis = paste0("_", version.Data, "_20210301")

## Directories to save results 
design.file = "../exp_design/RNAseq_sampleInfos.xlsx"
#dataDir = "../data/quantseq/featurecounts_R7183"

resDir = paste0("../results/", version.Data)
tabDir =  paste0(resDir, "/tables/")
tfDir = '~/workspace/imp/positional_memory/results/motif_analysis'
RdataDir = paste0(resDir, "/Rdata/")

if(!dir.exists(resDir)){dir.create(resDir)}
if(!dir.exists(tabDir)){dir.create(tabDir)}
if(!dir.exists(RdataDir)){dir.create(RdataDir)}

##########################################
# gene annotation 
##########################################
load(file=paste0(RdataDir, 'Design_stats_readCounts_', version.analysis, '.Rdata'))

annot = readRDS(paste0('/Volumes/groups/tanaka/People/current/jiwang/Genomes/axolotl/annotations/', 
                       'geneAnnotation_geneSymbols_cleaning_synteny_sameSymbols.hs.nr.rds'))

design$batch = as.factor(design$batch)

Select.genes.having.symbols = FALSE
if(Select.genes.having.symbols){
  gene.mapping = gene.mapping[which(!is.na(gene.mapping$gene.symbol.nr) | !is.na(gene.mapping$gene.symbol.hs)), ]
  all = all[!is.na(match(all$gene, gene.mapping$gene.id)), ]
}

QC.for.cpm = FALSE
if(QC.for.cpm){
  
  #index.qc = c(4)
  source(RNA.functions)
  source(RNA.QC.functions)
  
  kk = grep('Mature_UA|Mature_LA|Mature_Hand', design$conditions)
  kk = kk[which(design$batch[kk] == 1)]
  raw = all[, -1]
  raw = raw[, -kk]
  
  ss = apply(as.matrix(raw), 1, sum)
  
  
  raw = raw[which(ss >10), ]
  
  pdfname = paste0(resDir, "/Data_qulity_assessment_remove.MatureBatch1_filterGenes.withoutSymbols_10", version.analysis, ".pdf")
  pdf(pdfname, width = 12, height = 10)
  
  Check.RNAseq.Quality(read.count=raw, design.matrix = design[-kk, c(1, 2, 5)], 
                       lowlyExpressed.readCount.threshold=500)
  
  dev.off()
  
}

##########################################
# because the gene symbol from nr and hs are not consistent sometimes, so we keep gene.id from AMEXDD60
# Dimensionality reduction to visulize the difference between time points
# Here we select only the batch 3 and batch 2
##########################################
require(ggplot2)
require(DESeq2)
library("dplyr")
library("ggplot2")

raw = as.matrix(all[, -1])
rownames(raw) = all$gene

sels = which(design$batch != 1)

design.matrix = design[sels, ]
raw = raw[, sels]

dds <- DESeqDataSetFromMatrix(raw, DataFrame(design.matrix), design = ~ conditions)

ss = rowSums(counts(dds))

hist(log2(ss), breaks = 200, main = 'log2(sum of reads for each gene)')

cutoff.gene = 30
cat(length(which(ss > cutoff.gene)), 'genes selected \n')

dds <- dds[ss > cutoff.gene, ]

# normalization and dimensionality reduction
dds <- estimateSizeFactors(dds)
fpm = fpm(dds, robust = TRUE)

save(fpm, design, file = paste0(tfDir, '/RNAseq_fpm_fitered.cutoff.', cutoff.gene, '.Rdata'))

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('conditions', 'batch'), returnData = FALSE)
print(pca)

pca2save = as.data.frame(plotPCA(vsd, intgroup = c('conditions', 'batch'), returnData = TRUE))
ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= conditions, shape = batch))  + 
  geom_point(size=3) + 
  geom_text(hjust = 0.7, nudge_y = 1, size=2.5)

plot(ggp) + ggsave(paste0(resDir, "/PCAplot_batch2.batch3.pdf"), width=12, height = 8)

library("glmpca")
gpca <- glmpca(counts(dds), L=2)

gpca.dat <- gpca$factors
gpca.dat$conditions <- dds$conditions
gpca.dat$batch <- dds$batch
gpca.dat$name = rownames(gpca.dat)

ggplot(gpca.dat, aes(x = dim1, y = dim2, label = name, color = conditions, shape = batch)) +
  geom_point(size =3) + coord_fixed() + ggtitle("glmpca - Generalized PCA") +
  geom_text(hjust = 0.7, nudge_y = 1, size=2.5) + 
  ggsave(paste0(resDir, "/GPCAplot_batch2.batch3.pdf"), width=12, height = 8)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( vsd$conditions, vsd$batch, sep = " - " )
colnames(sampleDistMatrix) <- NULL

mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
mds$name = rownames(mds)

ggplot(mds, aes(x = `1`, y = `2`, label = name, color = conditions, shape = batch)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data") +
  geom_text(hjust = 0.7, nudge_y = 1, size=2.5) 

##########################################
# dynamic genes  
##########################################
dds <- estimateDispersions(dds)
plotDispEsts(dds, ymin = 10^-4)

dds <- nbinomLRT(dds, reduced = ~1 )
res <- results(dds)
