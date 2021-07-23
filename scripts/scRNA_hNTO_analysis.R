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
cells.sel = grep('No_Stretch_No_RASAG_Day_3|No_Stretch_RASAG_Day_5', design$orig.ident)
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
plotHighestExprs(sce, n=50) + fontsize

ave.counts <- calculateAverage(as.matrix(counts(sce)))
hist(log10(ave.counts), breaks=100, main="", col="grey80",
     xlab=expression(log[10]~"average count"))

num.cells <- nexprs(sce, byrow=TRUE)
hist(log10(num.cells), breaks=100, main="", col="grey80",
     xlab=expression(log[10]~"number cells"))

smoothScatter(log10(ave.counts), num.cells, ylab="Number of cells",
              xlab=expression(Log[10]~"average count"))

genes.to.keep <- num.cells > 20 & ave.counts >= 10^-2   # detected in >= 5 cells, ave.counts >=5 but not too high
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

sce = readRDS(file=paste0(RdataDir, 'scRNA_rawReadCounts_metadata_day3.day5.no.stretch_QCed_geneFiltered_scranNorm',
                          version.analysis, '.rds'))

sce = sce[, which(sce$orig.ident == 'No_Stretch_RASAG_Day_5')]

Use.scTransform = FALSE
if(!Use.scTransform){
  nt = as.Seurat(sce, counts = 'counts', data = 'logcounts', assay = "RNA") # scran normalized data were kept in Seurat
  
  nfeatures = 1000
  nt <- FindVariableFeatures(nt, selection.method = "vst", nfeatures = nfeatures)
  
  top10 <- head(VariableFeatures(nt), 10) # Identify the 10 most highly variable genes
  
  plot1 <- VariableFeaturePlot(nt)
  plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  CombinePlots(plots = list(plot1, plot2)) # plot variable features with and without labels
  
  nt = ScaleData(nt, features = rownames(nt))
  
}else{
  library(sctransform)
  #if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("glmGamPoi")
  nt <- SCTransform(nt, variable.features.n =  2000, verbose = FALSE)
  
}


cell.cycle.regression = FALSE
if(cell.cycle.regression){
  s.genes <- cc.genes.updated.2019$s.genes
  g2m.genes <- cc.genes.updated.2019$g2m.genes
  
  nt <- CellCycleScoring(nt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
                             
  nt <- RunPCA(nt, features = c(s.genes, g2m.genes))
  DimPlot(nt, reduction = 'pca')
  
  nt <- ScaleData(nt, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(nt))
  nt <- RunPCA(nt, features = c(s.genes, g2m.genes))
  DimPlot(nt, reduction = 'pca')
  
  RidgePlot(nt, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)
  
}

# new normalization from Seurat
# tried regress out the pct_counts_Mt but works less well
#nt <- SCTransform(object = nt, variable.features.n = nfeatures) 
nt <- RunPCA(object = nt, features = VariableFeatures(nt), verbose = FALSE, weight.by.var = FALSE)
ElbowPlot(nt, ndims = 50)


nt <- FindNeighbors(object = nt, reduction = "pca", k.param = 20, dims = 1:20)
nt <- FindClusters(nt, resolution = 0.5, algorithm = 3)

nb.pcs = 20; n.neighbors = 30; min.dist = 0.2;
nt <- RunUMAP(object = nt, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)

DimPlot(nt, reduction = "umap", group.by = 'orig.ident') + ggtitle(paste0('day3 and day5'))

DimPlot(nt, reduction = "umap", group.by = 'seurat_clusters') + ggtitle('scran normalization')

FeaturePlot(nt, features = c('FOXA2', 'FOXA1', 'LEF1', 'ID1', 'SPRY4'))
VlnPlot(nt, features = c("FOXA2", 'OLIG2', 'FGF8', 'NKX2.2', 'BMP7'),  pt.size = 0.2, ncol = 4)






