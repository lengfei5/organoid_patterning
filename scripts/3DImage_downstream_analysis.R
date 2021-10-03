##########################################################################
##########################################################################
# Project: Organoid patterning project  
# Script purpose: After the 3D image were segmented and parameters were extreacted by Hannah
# this script is to concantinate all readout of image processing and to do statistics as downstream analysis
# the output of image processing are in output/ folder 
# The hierachy of folders is organized by:
# Experiment (e.g. 210217_hNTdrugs3_0310outputs) > 
# Conditions (many, e.g. BMP4_SBonlyb_Statistics) > 
# Surfaces_1_Statistics (for cyst) and Surfaces_2_Statistics (for floorplate or FoxA2 positive clusters)
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Sat Apr 10 16:22:40 2021
##########################################################################
##########################################################################

########################################################
########################################################
# Section I : concatinate the readouts of 3D image analysis
# map the floorplate with the cysts 
########################################################
########################################################
rm(list=ls())

# specific input and output folders
resDir = '../results/CellProfiler_LDN.timeSeries'
tabDir = paste0(resDir, '/tables')
Rdata = paste0('../results/Rdata')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(tabDir)) dir.create(tabDir)
if(!dir.exists(Rdata)) dir.create(Rdata)

save.table.each.condition = FALSE
CellProfiler = TRUE

Manally.extract.metadata = TRUE
DAPI.channel = 'C1'

metadataCorrection = FALSE
dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/patterning_organoid/Images_data/Teresa_LDN_titration/batch.processed_d5/'

analysis.verison = 'Teresa.LDN.timeSeries_20210917_d5'

##########################################
# find associated fp for cyst at each condition from the processing output of CellProfiler
##########################################
#make_mergedTables_fromSegementation_Imaris()
#make_mergedTable_fromSegementation_cellProfiler()

# clean image metadata
image = read.csv(file = paste0(dataDir, 'image.csv'))
image = image[, 
              grep('ExecutionTime|MD5|Width|PathName|Metadata|ModuleError|ImageSet_ImageSet|Series_|ProcessingStatus|Channel_|Height|Frame', 
                   colnames(image), invert = TRUE)]

image = data.frame(image, stringsAsFactors = FALSE)
# drop the absoute path of image sources ONLY if all imges were frorm the same folder
image = image[, grep('URL_', colnames(image), invert = TRUE)]

# give each image an unique name using the DAPI channel
image$name = gsub(paste0('_isotropic_', DAPI.channel, '.tif'), '', as.character(image$FileName_DNA))

Dummy.imageNumber = image$ImageNumber[which(image$name == 'DUMMY')]
image = image[which(image$ImageNumber !=  Dummy.imageNumber), ]

# extract condition from image name and make design matrix
if(Manally.extract.metadata){
  image$condition = sapply(image$name, function(x) {
    xx = unlist(strsplit(as.character(x), '_')); 
    xx = xx[c((length(xx)-2):length(xx))]; 
    paste0(xx, collapse = '.')} )
  
  image$condition = gsub('.d5', '', image$condition)
  image$condition = paste0(image$condition, '.d5')
  image$condition = gsub('LDn', 'LDN', image$condition)
  
  if(metadataCorrection){
    design = readRDS(file = paste0(Rdata, '/perturbation_design_hNTdrugs3_0310.rds'))
    image$condition = NA
    
    for(n in 1:nrow(image))
    {
      kk = grep(image$name[n], design$Original.Image.Name)
      if(length(kk) != 1) cat('Error')
      image$condition[n] = design$condition[kk] 
    }
  }
}

# cyst
cyst = read.csv(file = paste0(dataDir, 'organoid.csv'))
colsToKeep = c('ImageNumber', 'ObjectNumber', 
               "AreaShape_Volume",  "AreaShape_SurfaceArea", 
               "AreaShape_Center_X", "AreaShape_Center_Y", "AreaShape_Center_Z", 
               "AreaShape_EquivalentDiameter", "AreaShape_MajorAxisLength", "AreaShape_MinorAxisLength",
               "Children_foxa2cluster_Count", 
               'Intensity_MeanIntensity_FOXA2', 'Intensity_IntegratedIntensity_FOXA2',
               "Intensity_IntegratedIntensity_Olig2", "Intensity_MeanIntensity_Olig2"
)
kk = match(colsToKeep, colnames(cyst))
if(any(is.na(kk))){
  cat('columns missing: \n', paste0(colsToKeep[which(is.na(kk))], collapse = '\n'), '\n')
  kk = kk[which(!is.na(kk))]
  
}
cyst = cyst[, kk]
cyst = cyst[which(cyst$ImageNumber != Dummy.imageNumber), ]

# foxA2 clusters
fp = read.csv(file = paste0(dataDir, 'foxa2cluster.csv'))

if(length(which(fp$Children_RelateObjects_Count != 1)) > 0 ){
  cat('Warning : -- some foxA2 clusters do not have parents in the table \n')
  
  fp = fp[which(fp$Children_RelateObjects_Count == 1), ]
}

colsToKeep_cluster = c('ImageNumber', 'ObjectNumber', "Parent_organoid", 
                       "AreaShape_Volume",  "AreaShape_SurfaceArea", "AreaShape_Center_X", "AreaShape_Center_Y", 
                       "AreaShape_Center_Z", 
                       "AreaShape_EquivalentDiameter", "AreaShape_MajorAxisLength", "AreaShape_MinorAxisLength",
                       "Distance_Centroid_organoid",  
                       'Intensity_MeanIntensity_FOXA2', 'Intensity_IntegratedIntensity_FOXA2'
)
jj = match(colsToKeep_cluster, colnames(fp))

if(any(is.na(jj))){
  cat('columns missing: \n', paste0(colsToKeep_cluster[which(is.na(jj))], collapse = '\n'), '\n')
  jj = jj[which(!is.na(jj))]
}

fp = fp[,jj]

fp = fp[which(fp$ImageNumber != Dummy.imageNumber), ]


# merge image information and cyst first and then merge cyst and fp
source('Functions_3dImage.R')
res = merge_image.cyst.fp_fromCellProfiler(image, cyst, fp)

saveRDS(res, file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_', analysis.verison, '.rds'))

########################################################
########################################################
# Section II: extract relevant parameters and compare across conditions
# 
########################################################
########################################################
res = readRDS(file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_', analysis.verison, '.rds'))
res = data.frame(res, stringsAsFactors = FALSE)

DoubleCheck.CP.surfaceArea.volume = FALSE

##########################################
# filter cyst or/and floorplates using global parameters
# QC plots
##########################################
require(ggplot2)
require(grid)
require(gridExtra)
library(tidyr)
library(dplyr)

res$condition = as.factor(res$condition)
res$ID_cyst = res$ID
res$ID_fp = paste0(res$ImageNumber_fp, '_', res$ObjectNumber_fp)
res$ID_fp[is.na(res$ObjectNumber_fp)] = NA

res$sphericity_cyst = pi^(1/3)*(6*res$AreaShape_Volume_cyst)^(2/3) / res$AreaShape_SurfaceArea_cyst*4/3 # mysterious factor 4/3
res$sphericity_fp = pi^(1/3)*(6*res$AreaShape_Volume_fp)^(2/3) / res$AreaShape_SurfaceArea_fp*4/3 

## Old filtering for Imaris segemetation
#res = res[which(res$Distance.to.Image.Border.XY.Unit.µm_Distance.to.Image.Border.XY.Img1_cyst > 0), ]
#res = res[which(res$Overlapped.Volume.Ratio_fp > 0.5 | is.na(res$Overlapped.Volume.Ratio_fp)), ]
#res = res[which(res$Sphericity.Unit._Sphericity_cyst > 0.8), ]

##########################################
## cyst filtering
##########################################
cond.id = paste0(res$condition, '_', res$ID_cyst)
mm = match(unique(cond.id), cond.id)
xx = res[mm, ]
xx$volume.log10 = log10(xx$AreaShape_Volume_cyst)

if(DoubleCheck.CP.surfaceArea.volume){ doubleCheck.CP.surfaceArea.volume(xx);}

p0 = as_tibble(xx) %>% 
  group_by(condition) %>% tally() %>%
  ggplot(aes(x = condition, y = n, fill = condition)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none")  + 
  ggtitle('nb of cysts ') + 
  theme(axis.text.x = element_text(angle = 90))
  
p1 = ggplot(xx, aes(x = condition, y=AreaShape_Volume_cyst, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
  ggtitle('cyst volume') + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90))

ggplot(xx, aes(x = condition, y=sphericity_cyst, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
  ggtitle('cyst volume') + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90))

p2 = ggplot(xx, aes(x = volume.log10)) +
  geom_histogram(binwidth = 0.1)

p3 = ggplot(xx, aes(x = sphericity_cyst)) +
  geom_histogram(binwidth = 0.01)

p23 = ggplot(xx, aes(x = volume.log10, y = sphericity_cyst, color = condition)) +
  geom_point(size = 1) +
  geom_hline(yintercept=0.85, colour = "red") + geom_vline(xintercept = 4, colour = "red")

sels = which(xx$volume.log10 >=3.5 & xx$sphericity_cyst >=0.9)
xx = xx[sels, ]

p4 = as_tibble(xx) %>% 
  group_by(condition) %>% tally() %>%
  ggplot(aes(x = condition, y = n, fill = condition)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none")  + 
  ggtitle('nb of cysts ') +
  theme(axis.text.x = element_text(angle = 90))

p5 = ggplot(xx, aes(x = condition, y=AreaShape_Volume_cyst, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
  ggtitle('cyst volume') + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90))

ggplot(xx, aes(x = condition, y=sphericity_cyst, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
  ggtitle('cyst volume') + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90))

pdfname = paste0(resDir, '/QC_cystFiltering_before_after_', analysis.verison, '.pdf')
while (!is.null(dev.list()))  dev.off()
pdf(pdfname,  width = 20, height = 16)

grid.arrange(p0, p1, ncol = 1) 
grid.arrange(p2, p3, p23, nrow = 2) 
grid.arrange(p4, p5, ncol = 1) 

dev.off()

res = res[!is.na(match(res$ID_cyst, xx$ID_cyst)), ]

saveRDS(res, file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_cystFilering_', analysis.verison, '.rds'))

##########################################
# fp filtering
##########################################
res = readRDS(file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_cystFilering_', analysis.verison, '.rds'))

#params = extract.turing.parameters.cellProfiler(res)
cond.id = paste0(res$condition, '_', res$ID_fp)
mm = match(unique(cond.id), cond.id)
xx = res[mm, ]
xx = xx[!is.na(xx$ID_fp), ]

xx$volume.log10 = log10(xx$AreaShape_Volume_fp)
xx$volume.ratio = (xx$AreaShape_Volume_fp/xx$AreaShape_Volume_cyst)
xx$volume.ratio.log10 = log10(xx$AreaShape_Volume_fp/xx$AreaShape_Volume_cyst)
xx$dist.cyst.radius.ratio = xx$Distance_Centroid_organoid_fp/(xx$AreaShape_EquivalentDiameter_cyst/2)

ggplot(xx, aes(x = volume.log10, y = volume.ratio.log10)) +
  geom_point(size = 0.2) + 
  geom_hline(yintercept = log10(0.002), colour = "red") + 
  geom_vline(xintercept = 2, colour = "red") + 
  geom_vline(xintercept = 2.5, colour = "red")

ggplot(xx, aes(x = condition, y=volume.ratio.log10, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.1) + 
  ggtitle('volume ratio fp/cyst') + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, size = 10))

p1 = ggplot(xx, aes(x = condition, y=AreaShape_Volume_fp, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
  ggtitle('fp volume') + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, size = 10)) +
  scale_y_continuous(trans='log10')

p2 = ggplot(xx, aes(x = condition, y=sphericity_fp, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
  ggtitle('fp sphericity') + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, size = 10))

p3 = ggplot(xx, aes(x = condition, y=Intensity_MeanIntensity_FOXA2_fp, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.1) + 
  ggtitle('fp mean FoxA2 intensity') + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, size = 10))

p4 = ggplot(xx, aes(x = volume.log10)) +
  geom_histogram(binwidth = 0.1)

p5 = ggplot(xx, aes(x = sphericity_fp)) +
  geom_histogram(binwidth = 0.02)

p6 = ggplot(xx, aes(x = Intensity_MeanIntensity_FOXA2_fp)) +
  geom_histogram(binwidth = 0.001) +
  geom_vline(xintercept = 0.01, colour = "red")

p45 = ggplot(xx, aes(x = volume.log10, y = sphericity_fp)) +
  geom_point(size = 0.2) +
  #geom_hline(yintercepst=0.85, colour = "red") + 
  geom_vline(xintercept = 1.5, colour = "red") + 
  geom_vline(xintercept = 2, colour = "red")

p56 = ggplot(xx, aes(x = volume.log10, y = Intensity_MeanIntensity_FOXA2_fp)) +
  geom_point(size = 0.2) +
  #geom_hline(yintercept=0.01, colour = "red") + 
  geom_vline(xintercept = 1.5, colour = "red") + 
  geom_vline(xintercept = 2, colour = "red")

#xx = xx[which(xx$volume.log10 >=2 & xx$volume.ratio>=10^-2.5), ]

xx = xx[which(xx$volume.log10 > 2.5 & xx$Intensity_MeanIntensity_FOXA2_fp >0.01 &
                xx$Intensity_MeanIntensity_FOXA2_fp & 
                xx$volume.ratio <= 1.0 ), ]

p7 = ggplot(xx, aes(x = condition, y=volume.log10, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
  ggtitle('fp volume') + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, size = 10))

p8 = ggplot(xx, aes(x = condition, y=sphericity_fp, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
  ggtitle('fp sphericity') + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, size = 10))

p9 = ggplot(xx, aes(x = condition, y=Intensity_MeanIntensity_FOXA2_fp, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.1) + 
  ggtitle('fp mean FoxA2 intensity') + theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, size = 10))


pdfname = paste0(resDir, '/QC_plots_fp_volumeFiltering_beforeAndAfter', analysis.verison, '.pdf')
while (!is.null(dev.list()))  dev.off()
pdf(pdfname,  width = 20, height = 20)

grid.arrange(p0, p1, p3, ncol = 1) 
grid.arrange(p4, p5, p6, nrow = 2)
grid.arrange(p45, p56, nrow = 2)
grid.arrange(p7, p8, p9, ncol = 1) 

dev.off()

#res$AreaShape_Volume_fp.log10 = log10(res$AreaShape_Volume_fp)
res = res[which(is.na(res$ID_fp) | !is.na(match(res$ID_fp, xx$ID_fp))), ]

saveRDS(res, file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_cyst.fp.Filering_', analysis.verison, '.rds'))

##########################################
# extract turing-relevant parameters
######################################### 
#res = readRDS(file = '~/workspace/imp/organoid_patterning/results/Rdata/mergedTable_cyst.fp_allConditions_cyst.fp.Filering_hNTdrugs3_0310_Analysis20210507_cellProfiler.rds')
#res = readRDS(file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_cyst.fp.Filering_', analysis.verison, '.rds'))

source('Functions_3dImage.R')
cat(length(unique(res$ID_cyst)), ' cysts and ', length(unique(res$ID_fp)) -1, 'fps\n')

params = extract.turing.parameters.cellProfiler(res, pixel.scale = 3, cyst.overlapRatio.threshold = 0.01)

saveRDS(params, file = paste0(Rdata, '/turing_parameters_extracted_', analysis.verison, '.rds'))

########################################################
########################################################
# Section : visualize the turing-model relevant parameters
#  
########################################################
########################################################
# merge results from multiple cellprofiler runs
params = readRDS(file = paste0(Rdata, '/turing_parameters_extracted_Teresa.LDN.timeSeries_20210917_d6.control.rds'))
rownames(params) = paste0('a', rownames(params))
xx = readRDS(file = paste0(Rdata, '/turing_parameters_extracted_Teresa.LDN.timeSeries_20210917_d6.LDN.rds'))
rownames(xx) = paste0('b', rownames(xx))

params = rbind(params, xx)

xx = readRDS(file = paste0(Rdata, '/turing_parameters_extracted_Teresa.LDN.timeSeries_20210917_d5.rds'))
rownames(xx) = paste0('c', rownames(xx))
params = rbind(params, xx)

saveRDS(params, file = paste0(Rdata, '/turing_parameters_extracted_Teresa.control.LDN_d5.d6_merged.rds'))


# read the merged table
params = readRDS(file = paste0(Rdata, '/turing_parameters_extracted_Teresa.control.LDN_d5.d6_merged.rds'))
params$condition = gsub('.d6.d6', '.d6', params$condition)


#library(tidyquant)
conds = unique(params$condition)

conds.sels = list(
   which(params$condition == "RA_LDNSB"))
  
# check controls first
conds.sels = list(
  # which(params$condition == "RA_LDNSB"|
  #             params$condition == 'noRA_LDNSB'| 
  #             params$condition == 'RA_LDNonly'| # controls 
  #         params$condition == 'RA_noLDNnoSB'),
  
  # which(params$condition == "RA_LDNSB"| # control
  #         params$condition == 'FGF_LDNSB'| # FGF single perburbation
  #         params$condition == 'PD_LDNSB'|
  #         params$condition == 'BMP4_SBonlyb'| # BMP single perturbation
  #         params$condition == 'RA_SBonly'|
  #         params$condition == 'FGF_SBonlyb'| # BMP FGF perturbation
  #         params$condition == 'PD_SBonly'
  #         )
  grep('LDN', params$condition, invert = TRUE), 
  setdiff(grep('d6', params$condition), 
          grep('N2B27', params$condition))
            
  
)

saveTalbe = FALSE

if(saveTalbe){
  nb.fp = as.numeric(as.character(params$nb.fp))
  sels = which(nb.fp>=0 & nb.fp<=10)
  xx = as.data.frame(as_tibble(params[sels, ]) %>% 
    group_by(condition, nb.fp) %>% tally())
  yy = matrix(0, ncol=6,  nrow = length(unique(xx$condition)))
  rownames(yy) = unique(xx$condition)
  colnames(yy) = paste0('foxa2Cluster.nb.', c(0:(ncol(yy)-1)))
  
  for(n in 1:nrow(yy))
  {
    for(kk in 0:5)
    {
      ii = which(xx$condition == rownames(yy)[n] & xx$nb.fp == kk)
      if(length(ii) == 1) yy[n, kk] = xx$n[ii]
    }
  }
  yy = data.frame(yy)
  yy$total.nb.cyst = apply(as.matrix(yy), 1, sum)
  
  yy = yy[which(yy$total.nb.cyst > 5), ]
  
  write.csv(yy, file = paste0(resDir, '/distribution_foxaCluster.nb.per.cyst_acrossConditions.csv'))
  
}

pdfname = paste0(resDir, '/NTorganoid_mouse_Teresa.d5.d6.contro.LDN.titration.pdf')
pdf(pdfname,  width = 20, height = 16)

for(n in 1:length(conds.sels))
{
  # n = 1
  sels = conds.sels[[n]]
  nb.fp = as.numeric(as.character(params$nb.fp[sels]))
  sels = sels[which(nb.fp>=0 & nb.fp<6)]
  
  if(n == 1){
    level_order = levels(factor(params$condition[sels]))
  }
  if(n == 2){
    level_order = c("F4.normRA.d6", "R1.normRA.d6", "R1.advRA.d6",  "R1.advRA.LDN10nM.d6", 
                    "R1.advRA.LDN50nM.d6", "R1.advRA.LDN100nM.d6", "R1.advRA.LDN200nM.d6", 
                    "R1.advRA.LDN300nM.d6", "R1.advRA.LDN400nM.d6", "R1.advRA.LDN500nM.d6", 
                    "R1.advRA.LDN750nM.d6", "R1.advRA.LDN1000nM.d6")
  }
  
  p0 = as_tibble(params[sels, ]) %>% 
    group_by(condition, nb.fp) %>% tally() %>%
    ggplot(aes(x = factor(condition, levels = level_order), y = n, fill = nb.fp)) +
    geom_bar(stat = "identity") +
    theme_classic() + ggtitle('nb of cysts and fp nb distribution') +
    theme(axis.text.x = element_text(angle = 90, size = 12))
  
  p1 = ggplot(params[sels, ], aes(x = factor(condition, levels = level_order), y=volume, fill=condition)) + 
    geom_violin(width = 2.0) + ggtitle('cyst volume') +
    theme(axis.text.x = element_text(angle = 90, size = 12))
    
  p2 = ggplot(params[sels, ], aes(x = factor(condition, levels =level_order), y=overlap.ratio, fill=condition)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.1) +  ylim(0, 0.75) + 
    ggtitle('cyst fraction overlapped by fp') +
    theme(axis.text.x = element_text(angle = 90, size = 12))
  
  p3 = ggplot(params[sels[which(as.numeric(as.character(params$nb.fp[sels]))>=1)], ], aes(x=factor(condition, levels = level_order), y=volume.fp, fill=condition)) + 
    geom_violin(width = 1.5) + ggtitle('foxa2 volume') +
    theme(axis.text.x = element_text(angle = 90, size = 12)) 
  
  p4 = ggplot(params[sels[which(as.numeric(as.character(params$nb.fp[sels]))>=1)], ], aes(x = factor(condition, levels = level_order), y=foxa2.fp, fill=condition)) + 
    geom_violin(width = 1.5) + ggtitle('FoxA2 mean intensity') + ylim(0, 20000) + 
    theme(axis.text.x = element_text(angle = 90, size = 12))
  
  p5 = ggplot(params[sels, ], aes(fill=condition, y=olig2 , x = condition)) + 
    geom_violin() + ggtitle('Olig2 mean intensity') +
    theme(axis.text.x = element_text(angle = 90, size = 12))
  
  #factor(Var1, levels = c('pooled', 'negative',  'positive'))))
  
  p6 = ggplot(params[sels, ], aes(x=nb.fp, y=volume, color=condition, fill = condition)) +
    geom_boxplot() + ggtitle('size dependency of fp nb (cyst volume)') + 
    theme(axis.text.x = element_text(angle = 0, size = 16), 
          axis.text.y = element_text(size = 16) )
  
  p61 = ggplot(params[sels, ], aes(x=nb.fp, y=radius.cyst, color=condition, fill = condition)) +
    geom_violin() + ggtitle('size dependency of fp nb (cyst radius)') 
  
  p7 = ggplot(params[sels[which(as.numeric(as.character(params$nb.fp[sels]))>1)], ], 
              aes(x=volume, y=dist.fp, color=condition)) +
    geom_point(size = 2.5) + ggtitle('distance between fps (wavelength)') 
   
  p8 = ggplot(params[sels[which(as.numeric(as.character(params$nb.fp[sels]))>1)], ], 
         aes(x=factor(condition, levels = level_order), y=dist.fp, color=condition, fill = condition)) +
    #geom_violin() + 
    geom_boxplot() +
    #geom_jitter(width = 0.1, color = 'black', size = 1.0) +
    ggtitle('distance between fps (wavelength)') +
    theme(axis.text.x = element_text(angle = 90, size = 12))
  
  p9 = ggplot(params[which(params$condition == 'R1.normRA.d6'|params$condition == 'R1.advRA.d6' ), ], aes(x=nb.fp, y=volume, color=condition, fill = condition)) +
    geom_boxplot() + ggtitle('size dependency of fp nb (cyst volume)') 
  
  grid.arrange(p0, p1,  nrow = 2, ncol = 1)
  grid.arrange(p2, p5,  nrow = 2, ncol = 1)
  grid.arrange(p7, p8,  nrow = 2, ncol = 1)
  grid.arrange(p4, p3,  nrow = 2, ncol = 1)
  grid.arrange(p6, p61,  nrow = 2, ncol = 1)
  grid.arrange(p9, nrow = 1, ncol = 1)
  
}

dev.off()
