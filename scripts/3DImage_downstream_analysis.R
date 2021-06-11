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
dataDir = '../images/210525_CellProfiler/'
resDir = '../results/210525_CellProfiler'
tabDir = paste0(resDir, '/tables')
Rdata = paste0('../results/Rdata')
analysis.verison = 'hNTdrugs3_0310_Analysis20210507_cellProfiler'

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(tabDir)) dir.create(tabDir)
if(!dir.exists(dataDir)) dir.create(dataDir)
if(!dir.exists(Rdata)) dir.create(Rdata)

save.table.each.condition = FALSE
CellProfiler = TRUE

metadataCorrection = TRUE

##########################################
# find associated fp for cyst at each condition
##########################################
if(CellProfiler){
  #make_mergedTable_fromSegementation_cellProfiler()
  
  # clean image metadata
  image = read.csv(file = paste0(dataDir, 'image.csv'))
  image = image[, 
  grep('ExecutionTime|MD5|Width|PathName|Metadata|ModuleError|ImageSet_ImageSet|Series_|ProcessingStatus|Channel_|Height|Frame', 
                       colnames(image), invert = TRUE)]
  
  image = data.frame(image, stringsAsFactors = FALSE)
  # drop the absoute path of image sources ONLY if all imges were frorm the same folder
  image = image[, grep('URL_', colnames(image), invert = TRUE)]
  
  # give each image an unique name
  image$name = gsub('_isotropic_C4.tif', '', as.character(image$FileName_DNA))
  
  Dummy.imageNumber = image$ImageNumber[which(image$name == 'DUMMY')]
  image = image[which(image$ImageNumber !=  Dummy.imageNumber), ]
  
  # extract condition from image name
  image$conds = sapply(image$name, function(x) {xx = unlist(strsplit(as.character(x), '_')); 
  xx = xx[-c(1, length(xx)-1, length(xx))]; paste0(xx, collapse = '.')} )
  
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
  cyst = cyst[, kk]
  cyst = cyst[which(cyst$ImageNumber != Dummy.imageNumber), ]
  
  # foxA2 clusters
  fp = read.csv(file = paste0(dataDir, 'foxa2cluster.csv'))
  
  if(length(unique(fp$Children_RelateObjects_Count)) != 1  | unique(fp$Children_RelateObjects_Count) != 1){
    cat('Error : -- only foxA2 cluster with parents should be in the table \n')
  }
  
  colsToKeep_cluster = c('ImageNumber', 'ObjectNumber', "Parent_organoid", 
                 "AreaShape_Volume",  "AreaShape_SurfaceArea", "AreaShape_Center_X", "AreaShape_Center_Y", 
                 "AreaShape_Center_Z", 
                 "AreaShape_EquivalentDiameter", "AreaShape_MajorAxisLength", "AreaShape_MinorAxisLength",
                 "Distance_Centroid_organoid",  
                 'Intensity_MeanIntensity_FOXA2', 'Intensity_IntegratedIntensity_FOXA2'
  )
  
  jj = match(colsToKeep_cluster, colnames(fp))
  fp = fp[,jj]
  
  fp = fp[which(fp$ImageNumber != Dummy.imageNumber), ]
  
  if(metadataCorrection){
    design = readRDS(file = paste0(Rdata, '/perturbation_design_hNTdrugs3_0310.rds'))  
  }
  
  save(image, cyst, fp, file = paste0(Rdata, '/image_cyst_fp_', analysis.verison, '.Rdata'))
  
  ##########################################
  # merge image information and cyst first
  # and then merge cyst and fp
  ##########################################
  load(file = paste0(Rdata, '/image_cyst_fp_', analysis.verison, '.Rdata'))
  
  # change to data.frame
  cyst = data.frame(cyst, stringsAsFactors = FALSE)
  image = data.frame(image, stringsAsFactors = FALSE)
  fp = data.frame(fp, stringsAsFactors = FALSE)
  
  # add suffix for each table in the colnames
  colnames(image) = paste0(colnames(image), '.image')
  colnames(cyst) = paste0(colnames(cyst), '.cyst')
  colnames(fp) = paste0(colnames(fp), '.fp')
  
  # merge the image and cyst tables
  kk = match(cyst$ImageNumber.cyst, image$ImageNumber.image)
  cyst = data.frame(image[kk, ], cyst, stringsAsFactors = FALSE)
  
  # start to merge cyst and fp
  cyst$ID = paste0(cyst$ImageNumber.image, '_', cyst$ObjectNumber.cyst)
  fp$parentID = paste0(fp$ImageNumber.fp, '_', fp$Parent_organoid.fp)
  
  # this is part of whole table in which there are only cysts with fp
  jj = match(fp$parentID, cyst$ID)
  res = data.frame(cyst[jj, ], fp, stringsAsFactors = FALSE) 
  
  # this is the second part of whole table for cysts without fp
  kk = which(is.na(match(cyst$ID, fp$parentID)))
  xx = matrix(NA, nrow = length(kk), ncol = ncol(fp))
  colnames(xx) = colnames(fp)
  xx = data.frame(cyst[kk, ], xx, stringsAsFactors = FALSE)
  
  # combine two parts to have whole table
  res = data.frame(rbind(res, xx))
  
  # sort the table with image number and cyst number
  res = res[with(res, order(ImageNumber.cyst, ObjectNumber.cyst)),  ]
  
  saveRDS(res, file = paste0(Rdata, 'mergedTable_cyst.fp_allConditions_', analysis.verison, '.rds'))
  
}else{
  make_mergedTables_fromSegementation_Imaris()
  
}
########################################################
########################################################
# Section II: extract relevant parameters and compare across conditions
# 
########################################################
########################################################
res = readRDS(file = paste0(Rdata, 'mergedTable_cyst.fp_allConditions_', analysis.verison, '.rds'))
res = data.frame(res, stringsAsFactors = FALSE)

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
cond.id = paste0(res$condition, '_', res$ID_cyst)
mm = match(unique(cond.id), cond.id)
xx = res[mm, ]

p0 = as_tibble(res[mm, ]) %>% 
  group_by(condition) %>% tally() %>%
  ggplot(aes(x = condition, y = n, fill = condition)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none")  + 
  ggtitle('nb of cysts ')
aes(class, hwy)
p1 = ggplot(res[mm, ], aes(x = condition, y=Overlapped.Volume.Ratio_cyst, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
  ggtitle('cyst fraction overlapped by fp') + theme(legend.position = "none") 

p2 = ggplot(res[mm, ], aes(x = condition, y=Intensity.Mean.Unit._Intensity.Mean.Ch3.Img1_cyst, fill=condition)) + 
  geom_violin() + ggtitle('FoxA2 mean intensity') + theme(legend.position = "none")

pdfname = paste0(resDir, '/QC_plots_beforeFiltering.pdf')
pdf(pdfname,  width = 25, height = 18)
grid.arrange(p0, p1, p2, nrow = 3, ncol = 1) 
dev.off()

# filter the images touching the boarder
res = res[which(res$Distance.to.Image.Border.XY.Unit.µm_Distance.to.Image.Border.XY.Img1_cyst > 0), ]

cond.id = paste0(res$condition, '_', res$ID_cyst)
mm = match(unique(cond.id), cond.id)

px = ggplot(xx, aes(x = Distance.to.Image.Border.XY.Unit.µm_Distance.to.Image.Border.XY.Img1_cyst)) +
  geom_histogram(binwidth = 10)
pxx = ggplot(res[mm, ], aes(x = Distance.to.Image.Border.XY.Unit.µm_Distance.to.Image.Border.XY.Img1_cyst)) +
  geom_histogram(binwidth = 10)

xx = res[mm, ]

p0 = as_tibble(res[mm, ]) %>% 
  group_by(condition) %>% tally() %>%
  ggplot(aes(x = condition, y = n, fill = condition)) +
  geom_bar(stat = "identity") +
  theme(legend.position = "none")  + 
  ggtitle('nb of cysts ')

p1 = ggplot(res[mm, ], aes(x = condition, y=Overlapped.Volume.Ratio_cyst, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
  ggtitle('cyst fraction overlapped by fp') + theme(legend.position = "none") 

p2 = ggplot(res[mm, ], aes(x = condition, y=Intensity.Mean.Unit._Intensity.Mean.Ch3.Img1_cyst, fill=condition)) + 
  geom_violin() + ggtitle('FoxA2 mean intensity') + theme(legend.position = "none")


pdfname = paste0(resDir, '/QC_plots_filteringImageOnBoard.pdf')
pdf(pdfname,  width = 25, height = 25)
grid.arrange(px, pxx, p0, p1, p2, nrow = 5, ncol = 1) 
dev.off()


res = res[which(res$Overlapped.Volume.Ratio_fp > 0.5 | is.na(res$Overlapped.Volume.Ratio_fp)), ]
#res = res[which(res$Overlapped.Volume.Ratio_cyst > 0.05 | is.na(res$Overlapped.Volume.Ratio_cyst)), ]

res = res[which(res$Sphericity.Unit._Sphericity_cyst > 0.8), ]

##########################################
# extract turing-relevant parameters
#########################################    

source('orgnoid_functions.R')
params = extract.turing.parameters(res)

##########################################
# visualize the turing-model relevant parameters
##########################################

for(n in 1:ncol(params))
{
  if(colnames(params)[n] == 'condition'|colnames(params)[n] == 'nb.fp'){
    params[ ,n] = as.factor(params[,n])
  }else{
    params[,n] = as.numeric(params[,n])
  }
}

conds = unique(params$condition)

# check controls first
conds.sels = list(
  # controls 
  which(params$condition == "RA_LDNSB"|
              params$condition == 'noRA_LDNSB'| 
              params$condition == 'RA_LDNonly'|
          params$condition == 'RA_noLDNnoSB'),
  
  which(params$condition == "RA_LDNSB"| 
          params$condition == 'FGF_LDNSB'|
          params$condition == 'PD_LDNSB'), # Fgf perturbation 
  
  which(params$condition == "RA_LDNSB"|
          params$condition == 'Chir3_LDNSB'|
          params$condition == 'Chir6_LDNSB'|
          params$condition == 'IWP2_LDNSBb'|
          params$condition == 'XAV_LDNSB'), # Wnt perturbation
  
  which(params$condition == "RA_LDNSB"|
          params$condition == 'BMP4_SBonlyb'| # BMP perturbation
          params$condition == 'RA_SBonly'
          )
  
  # which(params$condition == "RA_LDNSB"|
  #         params$condition == 'XAV_SBonly'| # BMP perturbation
  #         params$condition == 'IWP2_SBonly'|
  #         params$condition == 'Chir3_SBonly'|
  #         params$condition == 'Chir6_SBonly'
  # ), 
  # 
  # which(params$condition == "RA_LDNSB"|
  #         params$condition == 'FGF_SBonlyb'| # BMP perturbation
  #         params$condition == 'PD_SBonly'|
  #         params$condition == 'Chir3_SBonly'|
  #         params$condition == 'Chir6_SBonly'
  # )
  
)


pdfname = paste0(resDir, '/Organoid_perturbation_singlePerturbation.pdf')
pdf(pdfname,  width = 20, height = 10)

for(n in 1:length(conds.sels))
{
  # n = 1
  sels = conds.sels[[n]]
  nb.fp = as.numeric(as.character(params$nb.fp[sels]))
  sels = sels[which(nb.fp>=0 & nb.fp<10)]
  
  xx = params[sels, ]
  xx = xx[which(as.numeric(as.character(xx$nb.fp))> 1), ]
  xx = xx[which(xx$dist.fp > 300), ]
  
  #sels = c(1:nrow(params))
  #nb.fp = as.numeric(as.character(params$nb.fp[sels]))
  #sels = sels[which(nb.fp>=0 & nb.fp<10)]
  
  # general overview
  # p0 = as_tibble(params[sels, ]) %>% group_by(condition) %>% tally() %>%
  #   ggplot(aes(x = condition, y = n, fill = condition)) +
  #   geom_bar(stat = "identity") +
  #   theme_classic() + ggtitle('nb of cyst')
  
  p0 = as_tibble(params[sels, ]) %>% 
    group_by(condition, nb.fp) %>% tally() %>%
    ggplot(aes(x = condition, y = n, fill = nb.fp)) +
    geom_bar(stat = "identity") +
    theme_classic() + ggtitle('nb of cysts and fp nb distribution')
  
  p1 = ggplot(params[sels, ], aes(x = condition, y=volume, fill=condition)) + 
    geom_violin() + ggtitle('cyst volume')
  
  p2 = ggplot(params[sels, ], aes(x = condition, y=overlap.ratio, fill=condition)) + 
    geom_violin() + ggtitle('cyst fraction overlapped by fp')
  
  p3 = ggplot(params[sels, ], aes(x=condition, y=radius.fp, fill=condition)) + 
    geom_violin() + ggtitle('foxa2 radius') 
  
  p4 = ggplot(params[sels, ], aes(x = condition, y=foxa2.fp, fill=condition)) + 
    geom_violin() + ggtitle('FoxA2 mean intensity')
  
  p5 = ggplot(params[sels, ], aes(fill=condition, y=olig2 , x = condition)) + 
    geom_violin() + ggtitle('Olig2 mean intensity')
  
  # parameters relevant to Turing model
  #ggplot(params[sels, ], aes(x=area, y = radius.cyst, color=condition, fill = condition)) +
  #  geom_point() + ggtitle('estimated cyst radius')
  
  #ggplot(params[sels, ], aes(x=area, y = nb.fp, color=condition, fill = condition)) +
  
  #p1 = ggplot(xx, aes(x=radius.cyst.group, y=nb.fp, fill=condition)) +
  #  geom_violin() + ggtitle('estimated cyst radius')
  
  p6 = ggplot(params[sels, ], aes(x=nb.fp, y=radius.cyst, color=condition, fill = condition)) +
    geom_violin() + ggtitle('estimated cyst radius')
  
  p7 = ggplot(params[sels[which(as.numeric(as.character(params$nb.fp[sels]))>1)], ], aes(x=volume, y=dist.fp, color=condition)) +
    geom_point(size = 2.5) + ggtitle('distance between fps (wavelength)')
  
  p8 = ggplot(params[sels[which(as.numeric(as.character(params$nb.fp[sels]))>1)], ], aes(fill=condition, y=dist.fp, x = condition)) +
    geom_violin() + ggtitle('distance between fps (wavelength)')
  
  grid.arrange(p0, p6, p7, p8, p1, p2, p3, p4, p5,  nrow = 3, ncol = 3)
  
  #grid.arrange(p2, p3, nrow = 1, ncol = 2)
  
}

dev.off()


pdfname = paste0(resDir, '/Organoid_perturbation_summary_allConditions_v3.pdf')
pdf(pdfname,  width = 24, height = 8)

sels = c(1:nrow(params))
nb.fp = as.numeric(as.character(params$nb.fp[sels]))
sels = sels[which(nb.fp>=0 & nb.fp<10)]

p0 = as_tibble(params[sels, ]) %>% 
  group_by(condition, nb.fp) %>% tally() %>%
  ggplot(aes(x = condition, y = n, fill = nb.fp)) +
  geom_bar(stat = "identity") +
  theme_classic() + ggtitle('nb of cysts and fp nb distribution')

p1 = ggplot(params[sels, ], aes(x = condition, y=volume, fill=condition)) + 
  geom_violin() + ggtitle('cyst volume') + theme(legend.position = "none") 

p2 = ggplot(params[sels, ], aes(x = condition, y=overlap.ratio, fill=condition)) + 
  geom_violin(width=2) + ggtitle('cyst fraction overlapped by fp') + theme(legend.position = "none") +
  coord_cartesian(ylim = c(0, 0.1))

p3 = ggplot(params[sels, ], aes(x=condition, y=radius.fp, fill=condition)) + 
  geom_violin() + ggtitle('foxa2 radius') + coord_cartesian(ylim = c(0, 2*10^5)) + theme(legend.position = "none")

p4 = ggplot(params[sels, ], aes(x = condition, y=foxa2.fp, fill=condition)) + 
  geom_violin() + ggtitle('FoxA2 mean intensity') + theme(legend.position = "none")

p5 = ggplot(params[sels, ], aes(fill=condition, y=olig2 , x = condition)) + 
  geom_violin() + ggtitle('Olig2 mean intensity') + theme(legend.position = "none")

plot(p0)
plot(p1)
plot(p2)
plot(p3)
plot(p4)
plot(p5)

dev.off()


