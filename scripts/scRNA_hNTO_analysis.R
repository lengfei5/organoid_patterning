##########################################################################
##########################################################################
# Project: organoid patterning project
# Script purpose: analyze the scRNA-seq data from Adrian Ranga https://www.nature.com/articles/s41467-021-22952-0
# Here the Day 3 (before RA and SAG) and Day 5 (right after RA washing and without RA and SAG) data were selected;
# Day 11 data is not being considered, becasue the pattern at day 11 is not clear as salt-and-pepper, because it does not seem to be 
# stable pattern. 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Jul 22 13:50:47 2021
##########################################################################
##########################################################################
rm(list = ls())
RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
scRNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/scRNA_functions.R'
source(RNA.functions)
source(RNA.QC.functions)
source(scRNA.functions)

# setup for data import and sequencing QCs
version.analysis = '_v20210722'

resDir = paste0("../results/scRNAseq_hNTO_Adrian", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '../data/scRNA_Ranga/'

########################################################
########################################################
# Section : import the count table and metadata
# QCs 
########################################################
########################################################
library(data.table)
design = data.frame(fread(paste0(dataDir, 'GSE154120_hNTO_metadata.txt'), header=TRUE, sep="\t", stringsAsFactors=FALSE), 
                    stringsAsFactors = FALSE)
rownames(design) = design[, 1]
design = design[, -1]

raw = fread(paste0(dataDir, 'GSE154120_hNTO_processed_data.txt'), 
            header='auto', sep="\t", stringsAsFactors=FALSE)
raw = data.frame(raw, stringsAsFactors = FALSE)
rownames(raw) = raw[, 1]
raw = raw[, -1]

kk = match(rownames(design), colnames(raw))
raw = raw[, kk]

# select no stretch day 3 and day5 but not day11
cells.sel = grep('No_Stretch_No_RASAG_Day_3|No_Stretch_RASAG_Day_5|No_Stretch_RASAG_Day_11', design$orig.ident)
design = design[cells.sel, ]
raw = raw[, cells.sel]


save(raw, design, file=paste0(RdataDir, 'scRNA_rawReadCounts_metadata_day3.day5.no.stretch.Rdata'))

##########################################
# Check QC and the imported data was already filtered by the author
##########################################
library(SingleCellExperiment)
library(scater)
library(scRNA.seq.funcs)
library(scran)
options(stringsAsFactors = FALSE)

load(file=paste0(RdataDir, 'scRNA_rawReadCounts_metadata_day3.day5.no.stretch.Rdata'))
ss = apply(as.matrix(raw), 1, sum)
length(which(ss == 0))

options(stringsAsFactors = FALSE)

#load(file=paste0(RdataDir, version.DATA, '_RAW_Read_Counts_design_technicalRepMerged.Rdata'))
sce <- SingleCellExperiment(assays = list(counts = as.matrix(raw)), # count table should be a matrix
                            colData = as.data.frame(design),
                            rowData = data.frame(gene_names = rownames(raw), feature_symbol = rownames(raw)))


keep_feature <- rowSums(counts(sce) > 0) > 0
sce <- sce[keep_feature, ]

# check cell quality and no cell fitering
plotColData(sce, x="nCount_RNA", y="nFeature_RNA") + ggtitle("nFeatures vs umi counts") 
plotColData(sce, x="nCount_RNA", y="percent.mt") + ggtitle("% of Mt")

# plotColData(sce, y="log10_total_counts", x="seqInfos") + ggtitle("total nb of reads mapped to transcripts")
# plotColData(sce, y="total_features_by_counts", x="seqInfos") + ggtitle("total nb of genes")
# plotColData(sce,
#             x = "log10_total_counts",
#             y = "log10_total_features_by_counts",
#             #colour_by = "percent_mapped",
#             colour_by = "seqInfos",
#             size_by = "pct_counts_Mt"
# ) + scale_x_continuous(limits=c(4, 7)) +
#   scale_y_continuous(limits = c(2.5, 4.1)) +
#   geom_hline(yintercept=log10(c(500, 1000, 5000)) , linetype="dashed", color = "darkgray", size=0.5) +
#   geom_vline(xintercept = c(4:6), linetype="dotted", color = "black", size=0.5)


# filter genes
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotHighestExprs(sce, n=20) + fontsize

ave.counts <- calculateAverage(as.matrix(counts(sce)))
hist(log10(ave.counts), breaks=100, main="", col="grey80",
     xlab=expression(log[10]~"average count"))

num.cells <- nexprs(sce, byrow=TRUE)
hist(log10(num.cells), breaks=100, main="", col="grey80",
     xlab=expression(log[10]~"number cells"))

smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
              xlab=expression(Log[10]~"average count"))

genes.to.keep <- num.cells > 20 & ave.counts >= 10^-3   # detected in >= 5 cells, ave.counts >=5 but not too high
summary(genes.to.keep)

sce <- sce[genes.to.keep, ]

saveRDS(sce, file=paste0(RdataDir, 'scRNA_rawReadCounts_metadata_day3.day5.no.stretch_QCed_geneFiltered', version.analysis, '.rds'))


##########################################
# scran normalization 
##########################################
reducedDim(sce) <- NULL

sce$library.size = apply(counts(sce), 2, sum)
qclust <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = qclust)
sce <- logNormCounts(sce, log = TRUE, pseudo_count = 1)

par(mfrow = c(1, 1))
plot(sce$library.size/1e6, sizeFactors(sce), log="xy", xlab="Library size (millions)", ylab="Size factor", cex = 0.5)

saveRDS(sce, file=paste0(RdataDir, 'scRNA_rawReadCounts_metadata_day3.day5.no.stretch_QCed_geneFiltered_scranNorm',
                         version.analysis, '.rds'))


########################################################
########################################################
# Section : Seurat
# 
########################################################
########################################################
library(Seurat)
library(ggplot2)
library(dplyr)
library(SeuratWrappers)
library(cowplot)

sce = readRDS(file=paste0(RdataDir, 'scRNA_rawReadCounts_metadata_day3.day5.no.stretch_QCed_geneFiltered_scranNorm',
                          version.analysis, '.rds'))

#sce = sce[, which(sce$orig.ident == 'No_Stretch_RASAG_Day_5')]

Use.scTransform = FALSE
if(!Use.scTransform){
  nt = as.Seurat(sce, counts = 'counts', data = 'logcounts', assay = "RNA") # scran normalized data were kept in Seurat
  
  nfeatures = 5000
  nt <- FindVariableFeatures(nt, selection.method = "vst", nfeatures = nfeatures)
  nt = ScaleData(nt, features = rownames(nt))
  #top10 <- head(VariableFeatures(nt), 10) # Identify the 10 most highly variable genes
  #plot1 <- VariableFeaturePlot(nt)
  #plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  #CombinePlots(plots = list(plot1, plot2)) # plot variable features with and without labels
  
}else{
  library(sctransform)
  #if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("glmGamPoi")
  nt <- SCTransform(nt, variable.features.n =  2000, verbose = FALSE)
  
}

cell.cycle.regression = FALSE
if(cell.cycle.regression){
  ## see original code from https://satijalab.org/seurat/articles/cell_cycle_vignette.html
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  nt <- CellCycleScoring(nt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
                             
  nt <- RunPCA(nt, features = c(s.genes, g2m.genes))
  DimPlot(nt, reduction = 'pca')
  
  Regress.out.S.G2M = FALSE
  if(Regress.out.S.G2M){
    nt <- ScaleData(nt, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(nt))
    nt <- RunPCA(nt, features = c(s.genes, g2m.genes), verbose = FALSE)
    DimPlot(nt, reduction = 'pca')
    
    RidgePlot(nt, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
  }else{
    # As an alternative, we suggest regressing out the difference between the G2M and S phase scores
    nt$CC.Difference <- nt$S.Score - nt$G2M.Score
    nt <- ScaleData(nt, vars.to.regress = "CC.Difference", features = rownames(nt))
    
    nt <- RunPCA(nt, features = c(s.genes, g2m.genes), verbose = FALSE)
    DimPlot(nt, reduction = 'pca')
    
  }
  
  saveRDS(nt, file=paste0(RdataDir, 'scRNA_rawReadCounts_metadata_day3.day5.no.stretch_QCed_geneFiltered_scranNorm_cellCycleScoring',
                              version.analysis, '.rds'))
  
}

##########################################
# Run PCA, clusters and umap visualization
##########################################
nt = readRDS(file=paste0(RdataDir, 'scRNA_rawReadCounts_metadata_day3.day5.no.stretch_QCed_geneFiltered_scranNorm_cellCycleScoring',
                         version.analysis, '.rds'))

# subset nt for day5
#nt = subset(nt, cells = colnames(nt)[which(nt$orig.ident == 'No_Stretch_RASAG_Day_5')])
# nt = subset(nt, cells = colnames(nt)[which(nt$annotated_clusters == 'NP-5'|nt$annotated_clusters == 'V-5')])

Scran.HVGs = FALSE

# HVG with Seurat
nfeatures = 5000
nt <- FindVariableFeatures(nt, selection.method = "vst", nfeatures = nfeatures)

if(Scran.HVGs){
  
  sce = Seurat::as.SingleCellExperiment(nt, assay = 'RNA')
  dec <- modelGeneVar(sce)
  
  plot(dec$mean, dec$total, xlab="Mean log-expression", ylab="Variance", log = '', cex = 0.5)
  curve(metadata(dec)$trend(x), col="blue", add=TRUE)
  
  # Get the top 2000 genes.
  top.hvgs <- getTopHVGs(dec, n=nfeatures)
  
  # Get the top 10% of genes.
  top.hvgs2 <- getTopHVGs(dec, prop=0.1)
  
  # Get all genes with positive biological components.
  top.hvgs3 <- getTopHVGs(dec, var.threshold=0)
  
  # Get all genes with FDR below 5%.
  top.hvgs4 <- getTopHVGs(dec, fdr.threshold=0.05)
  
  
  dec = modelGeneCV2(sce)
  dec = dec[order(dec$p.value), ]
  top.hvgs.cv2 = rownames(dec)[1:nfeatures]
  
}

varFeatures = VariableFeatures(nt)
varFeatures = top.hvgs
varFeatures = top.hvgs.cv2

nt <- RunPCA(object = nt, features = varFeatures, verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(nt, ndims = 50)


nt <- FindNeighbors(object = nt, reduction = "pca", k.param = 20, dims = 1:30)
nt <- FindClusters(nt, resolution = 0.7, algorithm = 3)

nb.pcs = 30; n.neighbors = 20; min.dist = 0.2;
nt <- RunUMAP(object = nt, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)

DimPlot(nt, reduction = "umap", group.by = 'seurat_clusters') + ggtitle(paste0('all')) 

DimPlot(nt, reduction = "umap", group.by = 'orig.ident') + ggtitle(paste0('all')) 
DimPlot(nt, reduction = "umap", group.by = 'annotated_clusters') + ggtitle(paste0('all')) 

p1 = FeaturePlot(nt, features = c('FOXA2', 'SHH', 'ARX', 'HOXB4','ZNF703', 'CYP26A1', 'STRA8', 'HHIP1', 'PTCH1', 'GLI1', 'GLI2', 'GLI3'))
p2 = DimPlot(nt, reduction = "umap", group.by = 'orig.ident') + ggtitle(paste0('all')) 

CombinePlots(plots = list(p1, p2, ), ncol = 1)


xx = nt@assays$RNA@scale.data

kk = match(c('POU5F1', 'PAX6'), rownames(xx))

kk = match(c('SHH', 'FOXA2'), rownames(xx))

plot(xx[kk[1], ], xx[kk[2], ], xlab = 'SHH', ylab = 'FOXA2', cex  = 0.1)
jj = grep('Day_3', nt$orig.ident)
points(xx[kk[1],jj], xx[kk[2], jj], col = 'blue')

jj = grep('Day_5', nt$orig.ident)
points(xx[kk[1],jj], xx[kk[2], jj], col = 'red')

jj = grep('Day_11', nt$orig.ident)
points(xx[kk[1],jj], xx[kk[2], jj], col = 'green')


FeaturePlot(nt, features = c('FOXA2', 'OLIG2', 'NKX2-2', 'NKX2-1', 'NKX6-1', 'PAX6', 'SOX10', 'TPAP2C'))
#  ggsave(paste0(resDir, '/UMAP_day3_day5.pdf'), width = 12, height = 8)

#nt.sub = subset(nt, cells = colnames(nt)[which(nt$orig.ident == 'No_Stretch_RASAG_Day_5')])
#DimPlot(nt.sub, reduction = "umap", group.by = 'orig.ident') + ggtitle(paste0('day5'))
p1 = DimPlot(nt, reduction = "umap", group.by = 'seurat_clusters') + ggtitle('clusters with scran normalization')
p2 = FeaturePlot(nt, features = c("FOXA2", 'PAX6', 'NKX2-2', 'OLIG2'), ncol = 2)
p3 = VlnPlot(nt, features = c("FOXA2", 'PAX6', 'NKX2-2', 'OLIG2'), ncol = 2)
p4 = DimPlot(nt, reduction = "umap", group.by = 'annotated_clusters')
CombinePlots(plots = list(p1, p2, p3, p4), ncol = 2) + 
  ggsave(filename = paste0(resDir, '/markerGenes_for_FPcluster.pdf'), width = 16, height = 10)


examples = c('Lef1', 'Wnt3', 'Wnt3a', 'Wnt4', 'Wnt5b', 'Wnt6', 'Wnt7a', 'Wnt7b', 'Wnt8a', 'Dkk1', 'Dkk2', 'Dkk3', 'Tcf15', 'Tcf19', 
             "Wnt1", 'Sost', 'Sfrp5', 'Lypd6')
examples = c(c('Spry4', 'Spry2', 'Etv4', 'Etv5'), c('Fgf10', 'Fgf17','Fgf8', 'Fgf5', 'Fgf2','Fgf21', 'Fgf11','Fgf1', 'Fgf4', 'Fgfbp3'),
             'Dusp1', 'Dusp10', 'Dusp27', 'Dusp4', 'Dusp5')
examples = c('Id1', 'Id3', 'Smad6', 'Nog', 'Fst', 'Bambi', 'Bmp7', 'Bmp4', 'Bmp1', 'Bmp6', 'Bmpr2', 'Bmpr1b')

FeaturePlot(nt, features = c('FOXA2', toupper(examples)), ncol = 4)


nt$cells = 'Foxa2.neg'
nt$cells[which(nt$seurat_clusters == '5')] = 'Foxa2.pos'
Idents(nt) = 'cells'

saveRDS(nt, file = paste0(RdataDir, 
                          'scRNA_rawReadCounts_metadata_day5.no.stretch_QCed_geneFiltered_scranNorm_cellCycleScoring_FPclusterAnnoated',
                          version.analysis, '.rds'))


Use.subset.nt.to.annotate.FPcluster = FALSE
if(Use.subset.nt.to.annotate.FPcluster){
  #nt.sub = nt
  #saveRDS(nt.sub, file=paste0(RdataDir, 'scRNA_rawReadCounts_metadata_day5.no.stretch_QCed_geneFiltered_scranNorm_cellCycleScoring_NP5_V5',
  #                            version.analysis, '.rds'))
  nt.sub = readRDS(file=paste0(RdataDir, 'scRNA_rawReadCounts_metadata_day5.no.stretch_QCed_geneFiltered_scranNorm_cellCycleScoring_NP5_V5',
                               version.analysis, '.rds'))
  #nt = readRDS(file=paste0(RdataDir, 'scRNA_rawReadCounts_metadata_day3.day5.no.stretch_QCed_geneFiltered_scranNorm_cellCycleScoring',
  #                         version.analysis, '.rds'))
  #nt = subset(nt, cells = colnames(nt)[which(nt$orig.ident == 'No_Stretch_RASAG_Day_5')])
  
  nt$cells = 'Foxa2.neg'
  cell.pos = colnames(nt.sub)[which(nt.sub$cells == 'Foxa2.pos')]
  nt$cells[!is.na(match(colnames(nt), cell.pos))] = 'Foxa2.pos'
  Idents(nt) = 'cells'
  
}

nt = readRDS(file = paste0(RdataDir, 
                           'scRNA_rawReadCounts_metadata_day5.no.stretch_QCed_geneFiltered_scranNorm_cellCycleScoring_FPclusterAnnoated',
                           version.analysis, '.rds'))
table(nt$annotated_clusters)

nt = subset(nt, cells = colnames(nt)[which(nt$annotated_clusters != 'NC')])

DimPlot(nt, reduction = "umap", group.by = 'cells') + ggtitle(paste0('day5 Foxa2 positive and negative')) 

#VlnPlot(nt, features = c("FOXA2", 'OLIG2', 'FGF8', 'NKX2.2', 'BMP7', 'FGF4', 'SPRY4', 'ID1', 'LEF1', 
#                         'BMP4', 'FGF8', 'FGF2', 'FGF4', 'FGF10', 'WNT4'),  pt.size = 0.2, ncol = 4)

examples = c('Lef1', 'Wnt3', 'Wnt3a', 'Wnt4', 'Wnt5b', 'Wnt6', 'Wnt7a', 'Wnt7b', 'Wnt8a', 'Dkk1', 'Dkk2', 'Dkk3', 'Tcf15', 'Tcf19', 
             "Wnt1", 'Sost', 'Sfrp5', 'Lypd6')
VlnPlot(nt, features = c(toupper(examples)),  pt.size = 0.2, ncol = 4) + 
  ggsave(filename = paste0(resDir, '/FPcluster_otherClusters_WNT.pdf'), width = 16, height = 10)


examples = c('Spry4', 'Spry2', 'Etv4', 'Etv5', 'Fgf10', 'Fgf17','Fgf8', 'Fgf5', 'Fgf2','Fgf21', 'Fgf11','Fgf1', 'Fgf4', 'Fgfbp3',
             'Dusp1', 'Dusp10', 'Dusp27', 'Dusp4', 'Dusp5')
VlnPlot(nt, features = c(toupper(examples)),  pt.size = 0.2, ncol = 4)+ 
  ggsave(filename = paste0(resDir, '/FPcluster_otherClusters_FGF.pdf'), width = 16, height = 10)

examples = c('Id1', 'Id3', 'Smad6', 'Nog', 'Fst', 'Bambi', 'Bmp7', 'Bmp4', 'Bmp1', 'Bmp6', 'Bmpr2', 'Bmpr1b')
VlnPlot(nt, features = c(toupper(examples)),  pt.size = 0.2, ncol = 4)+ 
  ggsave(filename = paste0(resDir, '/FPcluster_otherClusters_BMP.pdf'), width = 16, height = 10)



##########################################
# aggreate the cell counts to have pseudo-bulk 
##########################################
sce = Seurat::as.SingleCellExperiment(nt, assay = 'RNA')

library(Matrix.utils)
# Subset metadata to only include the cluster and sample IDs to aggregate across
groups <- colData(sce)[, c("cells")]

# Aggregate across cluster-sample groups
out <- aggregate.Matrix(t(counts(sce)), groupings = groups, fun = "sum") 

out = t(out)

require(DESeq2)
require(ggplot2)
library(ggrepel)
require(gridExtra)
library(dplyr)
library(patchwork)
require(pheatmap)

raw = data.frame(out)
conds = colnames(raw)
dds <- DESeqDataSetFromMatrix(raw, DataFrame(conds = conds), design = ~ conds)

ss = rowSums(counts(dds))

hist(log2(ss), breaks = 200, main = 'log2(sum of reads for each gene)')

length(which(ss > quantile(ss, probs = 0.6)))
dd0 = dds[ss > quantile(ss, probs = 0.6) , ]
dd0 = estimateSizeFactors(dd0)
sizefactors.UQ = sizeFactors(dd0)

plot(sizeFactors(dd0), colSums(counts(dds)), log = 'xy')
text(sizeFactors(dd0), colSums(counts(dds)), colnames(dd0), cex =0.4)

sizeFactors(dds) = sizefactors.UQ
fpm = fpm(dds, robust = TRUE)

yy = data.frame(log2(fpm + 2^-6), gene = rownames(fpm), stringsAsFactors = FALSE)
yy$lfc = yy$Foxa2.pos - yy$Foxa2.neg

examples = c('Lef1', 'Wnt3', 'Wnt3a', 'Wnt4', 'Wnt5b', 'Wnt6', 'Wnt7a', 'Wnt7b', 'Wnt8a', 'Dkk1', 'Dkk2', 'Dkk3', 'Tcf15', 'Tcf19', 
             "Wnt1", 'Sost', 'Sfrp5', 'Lypd6')
examples = c(c('Spry4', 'Spry2', 'Etv4', 'Etv5'), c('Fgf10', 'Fgf17','Fgf8', 'Fgf5', 'Fgf2','Fgf21', 'Fgf11','Fgf1', 'Fgf4', 'Fgfbp3'),
             'Dusp1', 'Dusp10', 'Dusp27', 'Dusp4', 'Dusp5')
examples = c('Id1', 'Id3', 'Smad6', 'Nog', 'Fst', 'Bambi', 'Bmp7', 'Bmp4', 'Bmp1', 'Bmp6', 'Bmpr2', 'Bmpr1b')

examples = toupper(examples)
kk = c(); for(g in examples) kk = unique(c(kk, which(as.character(yy$gene) == g)))

ggplot(yy[kk, ], aes(x = gene, y = lfc)) + 
  geom_bar(stat = "identity", position="dodge") +
  coord_flip() + ggsave(paste0(resDir, "/LFC_positive.vs.negative_Day3.Day5_Wnt.pdf"), width=12, height = 10)

examples = c(c('Spry4', 'Spry2', 'Etv4', 'Etv5'), c('Fgf10', 'Fgf17','Fgf8', 'Fgf5', 'Fgf2','Fgf21', 'Fgf11','Fgf1', 'Fgf4', 'Fgfbp3'),
             'Dusp1', 'Dusp10', 'Dusp27', 'Dusp4', 'Dusp5')
kk = c(); for(g in examples) kk = unique(c(kk, which(as.character(yy$gene) == g)))
ggplot(yy[kk, ], aes(x = gene, y = lfc, fill = day)) + 
  geom_bar(stat = "identity", position="dodge") +
  coord_flip() + ggsave(paste0(resDir, "/LFC_positive.vs.negative_Day3.Day5_FGF.pdf"), width=12, height = 10)

examples = c('Id1', 'Id3', 'Smad6', 'Nog', 'Fst', 'Bambi', 'Bmp7', 'Bmp4', 'Bmp1', 'Bmp6', 'Bmpr2', 'Bmpr1b')

kk = c(); for(g in examples) kk = unique(c(kk, which(as.character(yy$gene) == g)))
ggplot(yy[kk, ], aes(x = gene, y = lfc, fill = day)) + 
  geom_bar(stat = "identity", position="dodge") +
  coord_flip() + ggsave(paste0(resDir, "/LFC_positive.vs.negative_Day3.Day5_BMP.pdf"), width=12, height = 10)


