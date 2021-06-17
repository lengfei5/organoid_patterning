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
  colnames(image) = paste0(colnames(image), '_image')
  colnames(cyst) = paste0(colnames(cyst), '_cyst')
  colnames(fp) = paste0(colnames(fp), '_fp')
  
  # merge the image and cyst tables
  kk = match(cyst$ImageNumber_cyst, image$ImageNumber_image)
  cyst = data.frame(image[kk, ], cyst, stringsAsFactors = FALSE)
  
  # start to merge cyst and fp
  cyst$ID = paste0(cyst$ImageNumber_image, '_', cyst$ObjectNumber_cyst)
  fp$parentID = paste0(fp$ImageNumber_fp, '_', fp$Parent_organoid_fp)
  
  index_cyst = c()
  index_fp = c()
  for(n in 1:nrow(cyst))
  {
    cat(n, '\n')
    kk = which(fp$parentID == cyst$ID[n])
    if(length(kk) == 0){
      index_cyst = c(index_cyst, n)
      index_fp = c(index_fp, NA)
    }else{
      index_cyst = c(index_cyst, rep(n, length(kk)))
      index_fp = c(index_fp, kk)
    }
  }
  
  res = data.frame(cyst[index_cyst, ], fp[index_fp, ], stringsAsFactors = FALSE) 
  # # this is part of whole table in which there are only cysts with fp
  # jj = match(fp$parentID, cyst$ID)
  # res = data.frame(cyst[jj, ], fp, stringsAsFactors = FALSE) 
  # 
  # # this is the second part of whole table for cysts without fp
  # kk = which(is.na(match(cyst$ID, fp$parentID)))
  # xx = matrix(NA, nrow = length(kk), ncol = ncol(fp))
  # colnames(xx) = colnames(fp)
  # xx = data.frame(cyst[kk, ], xx, stringsAsFactors = FALSE)
  # 
  # # combine two parts to have whole table
  # res = data.frame(rbind(res, xx))
  # 
  # # sort the table with image number and cyst number
  # res = res[with(res, order(ImageNumber_cyst, ObjectNumber_cyst)),  ]
  
  saveRDS(res, file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_', analysis.verison, '.rds'))
  
}else{
  make_mergedTables_fromSegementation_Imaris()
  
}
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

##########################################
## cyst filtering
##########################################
cond.id = paste0(res$condition, '_', res$ID_cyst)
mm = match(unique(cond.id), cond.id)
xx = res[mm, ]
xx$volume.log10 = log10(xx$AreaShape_Volume_cyst)

if(DoubleCheck.CP.surfaceArea.volume){
  #plot(xx$AreaShape_Volume_cyst, xx$AreaShape_SurfaceArea_cyst)
  #rr = c(0:1000)
  #points(4/3*pi*rr^3, 4*pi*rr^2, type = 'l')
  #points(rr^3, 6*rr^2, type = 'l')
  
  #plot(4*pi*(xx$AreaShape_EquivalentDiameter_cyst/2)^2, xx$AreaShape_SurfaceArea_cyst, cex = 0.6)
  #abline(0, 1, col = 'red')
  
  #plot(4/3*pi*(xx$AreaShape_EquivalentDiameter_cyst/2)^3, xx$AreaShape_Volume_cyst);
  #abline(0, 1, lwd =2.0, col = 'red')
}

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

p2 = ggplot(xx, aes(x = volume.log10)) +
  geom_histogram(binwidth = 0.1)

p3 = ggplot(xx, aes(x = sphericity_cyst)) +
  geom_histogram(binwidth = 0.01)

p23 = ggplot(xx, aes(x = volume.log10, y = sphericity_cyst)) +
  geom_point(size = 1) +
  geom_hline(yintercept=0.85, colour = "red") + geom_vline(xintercept = 4, colour = "red")

sels = which(xx$volume.log10 >=4 & xx$sphericity_cyst >=0.8)
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


pdfname = paste0(resDir, '/QC_cystFiltering_before_after.pdf')
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

p5 = ggplot(xx, aes(x = sphericity_cyst)) +
  geom_histogram(binwidth = 0.01)

p6 = ggplot(xx, aes(x = Intensity_MeanIntensity_FOXA2_fp)) +
  geom_histogram(binwidth = 0.001)

p45 = ggplot(xx, aes(x = volume.log10, y = sphericity_fp)) +
  geom_point(size = 0.2) +
  #geom_hline(yintercepst=0.85, colour = "red") + 
  geom_vline(xintercept = 1.5, colour = "red")

p56 = ggplot(xx, aes(x = volume.log10, y = Intensity_MeanIntensity_FOXA2_fp)) +
  geom_point(size = 0.2) +
  #geom_hline(yintercept=0.01, colour = "red") + 
  geom_vline(xintercept = 1.5, colour = "red")


sels = which(xx$volume.log10 >=1.5) 
xx = xx[sels, ]

p7 = ggplot(xx, aes(x = condition, y=AreaShape_Volume_fp, fill=condition)) + 
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


pdfname = paste0(resDir, '/QC_plots_fp_volumeFiltering_beforeAndAfter.pdf')
pdf(pdfname,  width = 20, height = 20)

grid.arrange(p0, p1, p3, ncol = 1) 
grid.arrange(p4, p5, p6, nrow = 2)
grid.arrange(p45, p56, nrow = 2)
grid.arrange(p7, p8, p9, ncol = 1) 

dev.off()

res$AreaShape_Volume_fp.log10 = log10(res$AreaShape_Volume_fp)
res = res[is.na(res$AreaShape_Volume_fp.log10) | res$AreaShape_Volume_fp.log10 >= 1.5, ]

saveRDS(res, file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_cyst.fp.Filering_', analysis.verison, '.rds'))


##########################################
# extract turing-relevant parameters
######################################### 
res = readRDS(file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_cyst.fp.Filering_', analysis.verison, '.rds'))
source('orgnoid_functions.R')

params = extract.turing.parameters.cellProfiler(res, pixel.scale = 3)

##########################################
# visualize the turing-model relevant parameters
##########################################
conds = unique(params$condition)

# check controls first
conds.sels = list(
  # controls 
  which(params$condition == "RA_LDNSB"|
              params$condition == 'noRA_LDNSB'| 
              params$condition == 'RA_LDNonly'|
          params$condition == 'RA_noLDNnoSB'),
  
  which(params$condition == "RA_LDNSB"| # control
          params$condition == 'FGF_LDNSB'| # FGF single perburbation
          params$condition == 'PD_LDNSB'|
          params$condition == 'BMP4_SBonlyb'| # BMP single perturbation
          params$condition == 'RA_SBonly'|
          params$condition == 'FGF_SBonlyb'| # BMP FGF perturbation
          params$condition == 'PD_SBonly'
          )
)


pdfname = paste0(resDir, '/Organoid_perturbation_singlePerturbation.pdf')
pdf(pdfname,  width = 20, height = 16)

for(n in 1:length(conds.sels))
{
  # n = 1
  sels = conds.sels[[n]]
  nb.fp = as.numeric(as.character(params$nb.fp[sels]))
  sels = sels[which(nb.fp>=0 & nb.fp<=10)]
  
  #xx = params[sels, ]
  #xx = xx[which(as.numeric(as.character(xx$nb.fp))> 1), ]
  #xx = xx[which(xx$dist.fp > 300), ]
  
  p0 = as_tibble(params[sels, ]) %>% 
    group_by(condition, nb.fp) %>% tally() %>%
    ggplot(aes(x = condition, y = n, fill = nb.fp)) +
    geom_bar(stat = "identity") +
    theme_classic() + ggtitle('nb of cysts and fp nb distribution') 
  
  p1 = ggplot(params[sels, ], aes(x = condition, y=volume, fill=condition)) + 
    geom_violin() + ggtitle('cyst volume') 
    
  p2 = ggplot(params[sels, ], aes(x = condition, y=overlap.ratio, fill=condition)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) + 
    ggtitle('cyst fraction overlapped by fp')
  
  p3 = ggplot(params[sels, ], aes(x=condition, y=radius.fp, fill=condition)) + 
    geom_violin() + ggtitle('foxa2 radius')
  
  p4 = ggplot(params[sels, ], aes(x = condition, y=foxa2.fp, fill=condition)) + 
    geom_violin() + ggtitle('FoxA2 mean intensity')
  
  p5 = ggplot(params[sels, ], aes(fill=condition, y=olig2 , x = condition)) + 
    geom_violin() + ggtitle('Olig2 mean intensity')
  
  p6 = ggplot(params[sels, ], aes(x=nb.fp, y=volume, color=condition, fill = condition)) +
    geom_violin() + ggtitle('size dependency of fp nb')
  
  p7 = ggplot(params[sels[which(as.numeric(as.character(params$nb.fp[sels]))>1)], ], aes(x=volume, y=dist.fp, color=condition)) +
    geom_point(size = 2.5) + ggtitle('distance between fps (wavelength)') 
   
  p8 = ggplot(params[sels[which(as.numeric(as.character(params$nb.fp[sels]))>1)], ], 
         aes(x=condition, y=dist.fp, color=condition, fill = condition)) +
    geom_violin() + ggtitle('distance between fps (wavelength)') 
  
    
  grid.arrange(p0, p1, p2, p5,  nrow = 2, ncol = 2)
  
  grid.arrange(p6, p7, p8, p4, p3,  nrow = 2, ncol = 3)
  
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


