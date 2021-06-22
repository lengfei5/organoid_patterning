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

Counts.to.Use = "UMI"
design = readRDS(file = paste0(RdataDir, 'sampleInfo_QC.stats.rds'))

colnames(design)[1] = 'SampleID'
design$conds = paste0(design$condition, '_', design$cells)

# prepare the data table 
Dir_umi = paste0(dataDir, "htseq_counts_BAMs_umi")
#Dir_read = paste0(dataDir, "htseq_counts_BAMs")

aa <- list.files(path = Dir_umi, pattern = "*umiDedup.txt", full.names = TRUE)
aa = merge.countTables.htseq(aa)
#colnames(aa1)[-1] = paste0(colnames(aa1)[-1], ".UMI")

#aa2 <- list.files(path = Dir_read, pattern = "*.txt", full.names = TRUE)
#aa2 = merge.countTables.htseq(aa2)
#colnames(aa2)[-1] = paste0(colnames(aa2)[-1], ".readCount")

#aa <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), list(aa1, aa2))

## compare read counts vs. umi counts
# source(RNAfunctions)
# Compare.UMI.vs.readCounts = TRUE
# if(Compare.UMI.vs.readCounts){
#   pdfname = paste0(resDir, "readCounts_vs_UMI_normalized", version.analysis, ".pdf")
#   pdf(pdfname, width = 10, height = 8)
#   
#   compare.readCount.UMI(design, aa, normalized = FALSE)
#   
#   dev.off()
# }

# if(Counts.to.Use == 'readCounts'){
#   all = process.countTable(all=aa, design = design, special.column = ".readCount", ensToGeneSymbol = TRUE)
# }else{
#   if(Counts.to.Use == "UMI"){
#     all = process.countTable(all=aa, design = design, special.column = "UMI", ensToGeneSymbol = TRUE)
#   }else{
#     cat("Error : no counts found for ", Counts.to.Use, "for miRNAs \n")
#   }
# }

all = process.countTable(all=aa, design = design[,c(1, 14)], special.column = NULL, ensToGeneSymbol = FALSE)

all = all[grep('^__', all$gene, invert = TRUE), ]

save(design, all, file=paste0(RdataDir, '/Design_Raw_readCounts', version.analysis, '.Rdata'))

########################################################
########################################################
# Section II : analyze the RNA-seq data 
# 
########################################################
########################################################

##########################################
# gene names converted from ensID to gene symbol 
##########################################
load(file=paste0(RdataDir, '/Design_Raw_readCounts', version.analysis, '.Rdata'))

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

save(design, all, file = paste0(RdataDir, '/design.detailed_RawUMI_', version.analysis, '.Rdata'))

##########################################
# QCs of replicates and conditions
##########################################
load(file = paste0(RdataDir, '/design.detailed_RawUMI_', version.analysis, '.Rdata'))

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

##########################################
# because the gene symbol from nr and hs are not consistent sometimes, so we keep gene.id from AMEXDD60
# Dimensionality reduction to visulize the difference between time points
# Here we select only the batch 3 and batch 2
##########################################
require(ggplot2)
require(DESeq2)
library("dplyr")

load(file = paste0(RdataDir, '/design.detailed_RawUMI_', version.analysis, '.Rdata'))

#raw = as.matrix(all[, -1])
#rownames(raw) = all$gene
#sels = which(design$batch != 1)
sels = which(design$SampleID != '161040' & design$SampleID != '161039')

design.matrix = design[sels, ]
raw = all[, sels]

dds <- DESeqDataSetFromMatrix(raw, DataFrame(design.matrix), design = ~ conds)

ss = rowSums(counts(dds))

hist(log2(ss), breaks = 200, main = 'log2(sum of reads for each gene)')

cutoff.gene = 100
cat(length(which(ss > cutoff.gene)), 'genes selected \n')

dds <- dds[ss > cutoff.gene, ]

# normalization and dimensionality reduction
dds <- estimateSizeFactors(dds)
fpm = fpm(dds, robust = TRUE)

#save(fpm, design, file = paste0(tfDir, '/RNAseq_fpm_fitered.cutoff.', cutoff.gene, '.Rdata'))

vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

pca=plotPCA(vsd, intgroup = c('condition', 'cells'), returnData = TRUE, ntop = 500)
#print(pca)
pca2save = as.data.frame(pca)

ggp = ggplot(data=pca2save, aes(PC1, PC2, label = name, color= condition, shape = cells))  + 
  geom_point(size=4) + 
  geom_text(hjust = 0.2, nudge_y = 0.4, size=3)

plot(ggp) + ggsave(paste0(resDir, "/PCAplot_ntop500.pdf"), width=12, height = 8)


Test.pooling.negative.positive.cells = FALSE
if(Test.pooling.negative.positive.cells){
  fpm = fpm(dds, robust = TRUE)
  design.matrix$condition.replicate = paste0(design.matrix$condition, '_', design.matrix$triplicate.nb)
  
  rrs = readRDS(file = paste0(RdataDir, 'facs_positive_negative_ratios.rds'))
  rrs$Tube.Name = as.character(rrs$Tube.Name)
  rrs$Tube.Name = gsub('[-]', '_', rrs$Tube.Name)
  rrs$Tube.Name = gsub('Fgf', 'FGF', rrs$Tube.Name)
  rrs$Tube.Name = gsub('PD03', 'PD', rrs$Tube.Name)
  rrs$Tube.Name = gsub('Chiron', 'CHIR', rrs$Tube.Name)
  
  pools = matrix(NA, nrow = nrow(fpm), ncol = length(unique(design.matrix$condition.rep)))
  colnames(pools) = unique(design.matrix$condition.rep)
  rownames(pools) = rownames(fpm)
  
  for(n in 1:ncol(pools))
  {
    # n = 1
    ratio = rrs$GFP_pos..Parent[which(rrs$Tube.Name == colnames(pools)[n])]/100  
    cat(n, ' - ', colnames(pools)[n], '- positive cell ratio ', ratio,  '\n')
    jj1 = which(design.matrix$condition.replicate == colnames(pools)[n] & design.matrix$cells == 'Foxa2.pos')
    jj2 = which(design.matrix$condition.replicate == colnames(pools)[n] & design.matrix$cells == 'Foxa2.neg')
    
    pools[,n] = ratio * fpm[,jj1] + (1 - ratio) * fpm[, jj2]
    
  }
  
  cc.pools = colnames(pools)
  cc.pools = sapply(cc.pools, function(x) unlist(strsplit(as.character(x), '_'))[1])
  
  library(factoextra)
  ntop = 3000
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
  
}
  
##########################################
# first check the normalized signals  
##########################################
#fpm = fpm(dds, robust = TRUE)
#fpm = log2(fpm + 2^-4)

fpm = pools

examples = unique(c('Foxa2', # FoxA 
                    'Lef1', 'Mapk1', rownames(fpm)[grep('Smad', rownames(fpm))], 
                    rownames(fpm)[grep('Wnt|Dkk|Tcf', rownames(fpm))],
                    rownames(fpm)[grep('Bmp', rownames(fpm))], 'Nog', 'Chrd', 'Runx1', 'Runx2',  'Smad6', 'Id1', 'Id3',
                    rownames(fpm)[grep('Acvr', rownames(fpm))],
                    rownames(fpm)[grep('Fgf', rownames(fpm))], 
                    'Dusp1', 'Dusp10', 'Dusp27', 'Dusp4', 'Dusp5', 'Mapk10', 'Mapk4', 'Mapk8ip2', 'Spry4', 'Rbpj', 'Hes1', 'Hes5',
                    'Hes7', 'Hey1', 'Hey2',
                    rownames(fpm)[grep('Notch|Jag|Dll|Dlk', rownames(fpm))], 
                    rownames(fpm)[grep('Tgf', rownames(fpm))]
))

#n = which(rownames(fpm) == 'Foxa2')
pdfname = paste0(resDir, '/TM3_examples_FoxA2.positive.vs.negative_POOLed_v1.pdf')
pdf(pdfname,  width = 10, height = 6)
par(cex = 1.0, las = 1, mgp = c(3,2,0), mar = c(6,6,2,0.2), tcl = -0.3)

reps = colnames(pools)
reps = sapply(reps, function(x) unlist(strsplit(as.character(x), '_'))[2])

for(g in examples)
{
  # g = 'Foxa2'
  kk = which(rownames(fpm) == g)
  
  if(length(kk) > 0){
    cat(g, '\n')
    
    #xx = data.frame(cpm = fpm[kk, ], condition = design.matrix$condition.rep,
    #                cells = design.matrix$cells, replicate= design.matrix$Triplicate.No.)
    
    xx = data.frame(cpm = fpm[kk, ], condition = cc.pools, rep = reps)
    
    p0 = ggplot(xx,  aes(x = condition, y = cpm, color = rep, fill = condition)) +
      geom_bar(stat = "identity", position="dodge") +
      geom_hline(yintercept = 2, colour = "blue") + 
      ggtitle(g)  + 
      theme(axis.text.x = element_text(angle = 90, size = 10))
    plot(p0)
    
    # plot(c(0, 1), type = 'n', xlim = c(-18, 60), ylim = range(c(fpm[kk, ], 0, 1)), main = g, 
    #      ylab = 'log2(fpm)', xlab = 'time')
    # points(tt, fpm.RA[kk, ], col = 'darkblue', type = 'l', pch = 16, lwd = 2.0)
    # points(tt, fpm.RA[kk, ], col = 'darkblue', type = 'p', pch = 16)
    # 
    # points(tt, fpm.noRA[kk, ], col = 'darkred', type = 'l', lwd = 1.0)
    # points(tt, fpm.noRA[kk, ], col = 'darkred', type = 'p', pch = 1)
    # 
    # points(48, sorted[kk, 1], col = 'darkblue', type = 'p', cex = 2.0, pch = 21, bg = 'magenta')
    # points(48, sorted[kk, 2], col = 'darkblue', type = 'p', cex = 2.0, pch = 21, bg = 'darkgreen')
    # #points(48, sorted[kk, 3], col = 'darkorange', type = 'p', cex = 2.0, pch = 18)
    # 
    #abline(h = c(0, 1), col = 'darkgray', lwd = 2.0)
    #legend('topleft', legend = c('RA', 'no.RA', 's48h.RA.AF', 's48h.RA.GFPp'), bty = 'n', 
    #       col = c('darkblue', 'darkred', 'magenta', 'darkgreen'), lwd =2.0, pch = c(16, 1, 16, 16), lty = c(1, 1, 0, 0))
    
  }else{
    cat(g, 'Not Found \n')
  }
}

dev.off()



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

