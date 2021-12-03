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

# load packages required (install missing packages)
require(grid)
require(gridExtra)
library(tidyr)
library(dplyr)
require(ggplot2)

source('Functions_3Dimage.R') # be careful of the directory of this function, better to put it in the same folder as script

# specific input and output folders
resDir = '../results/CellProfiler_Teresa_SAG.BMP4.titration'
analysis.verison = '20211130' # optional
tabDir = paste0(resDir, '/tables')
Rdata = paste0('../results/Rdata')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(tabDir)) dir.create(tabDir)
if(!dir.exists(Rdata)) dir.create(Rdata)

save.table.each.condition = FALSE
Manally.extract.metadata = TRUE
DAPI.channel = 'C1' # important to specify

dataDir = '/Volumes/groups/tanaka/People/current/Teresa/CellProfilerOutput/SAGBMP4titration_20211127/outputfiles_parentchild_20211127/'


##########################################
# find associated fp for cyst at each condition from the processing output of CellProfiler
# merge image information and cyst first and then merge cyst and fp
##########################################
#double check where you are, to print this info that you check in terminal
cat('input directory -- ', dataDir, '\n')

image = read.csv(file = paste0(dataDir, 'MyExpt_Image.csv'))

# call function 'clean_image_table'
image = clean_image_table(image, DAPI.channel)

# remove Dummy image if there is such image
cat('remove dummy image \n')
Dummy.imageNumber = NULL

if(length(which(image$name == 'DUMMY')) >0){
  if(!is.na(image$ImageNumber[which(image$name == 'DUMMY')])){
    Dummy.imageNumber = image$ImageNumber[which(image$name == 'DUMMY')]
    image = image[which(image$ImageNumber !=  Dummy.imageNumber), ]
  }
  
}else{
  cat('NO image wiht name DUMMY \n')
}

# cyst
cyst = read.csv(file = paste0(dataDir, 'MyExpt_organoid.csv'))
cyst = clean_cyst_table(cyst, Dummy.imageNumber)

# foxA2 clusters
fp = read.csv(file = paste0(dataDir, 'MyExpt_FOXA2cluster.csv'))
fp = clean_fp_table(fp, Dummy.imageNumber)


# extract condition from image name and make design matrix
# manually extract the condition, data-specific
# should do it for each for new dataset
if(Manally.extract.metadata){
  
  cat('specific condition for image \n')
  cat('should be specify for each new data \n')
  
  image$condition = sapply(image$name, function(x) {
    xx = unlist(strsplit(as.character(x), '_')); 
    #xx = xx[c((length(xx)-2):length(xx))]; 
    xx = xx[length(xx)]
    paste0(xx, collapse = '.')} )
  
  #image$condition = gsub('.d6', '', image$condition) # data all from d6
  #image$condition = paste0(image$condition, '.d6')
  #image$condition = gsub('LDn', 'LDN', image$condition)
  
}

res = merge_image.cyst.fp_fromCellProfiler(image, cyst, fp)

res = data.frame(res, stringsAsFactors = FALSE)

res = correct_sphericity_cyst.fp(res)
cat(range(res$sphericity_cyst, na.rm = TRUE))
cat(range(res$sphericity_fp, na.rm = TRUE))

saveRDS(res, file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_', analysis.verison, '.rds'))
rm(cyst, image, fp)

########################################################
########################################################
# Section II:  # filter cyst or/and floorplates using global parameters
# QC plots
# 
########################################################
########################################################

##########################################
## cyst filtering
##########################################
cond.id = paste0(res$condition, '_', res$ID_cyst)
cyst_all = res[match(unique(cond.id), cond.id), ]
cyst_all$volume.log10 = log10(cyst_all$AreaShape_Volume_cyst)


ggplot(cyst_all, aes(x = volume.log10, y = sphericity_cyst, color = condition)) +
  geom_point(size = 1) +
  geom_hline(yintercept=0.8, colour = "red") + geom_vline(xintercept = 3.5, colour = "red")

cat('Please choose the threshod for the cyst volume and sphericity for the cyst filtering \n')
threshold.cyst.volume = 3.5
threshold.cyst.sphericity = 0.8

sels = which(cyst_all$volume.log10 >= threshold.cyst.volume & cyst_all$sphericity_cyst >= threshold.cyst.sphericity)
cyst_sel = cyst_all[sels, ]

res = res[!is.na(match(res$ID_cyst, cyst_sel$ID_cyst)), ]

saveRDS(res, file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_cystFilering_', analysis.verison, '.rds'))

cat('Do you wannt to compare the parameters before and after cyst filtering')
Compare.cyst.fitering.before.vs.after = TRUE

if(Compare.cyst.fitering.before.vs.after){
  p0 = as_tibble(cyst_all) %>% 
    group_by(condition) %>% tally() %>%
    ggplot(aes(x = condition, y = n, fill = condition)) +
    geom_bar(stat = "identity") +
    theme(legend.position = "none")  + 
    ggtitle('nb of cysts ') + 
    theme(axis.text.x = element_text(angle = 90))
  
  p1 = ggplot(cyst_all, aes(x = condition, y=AreaShape_Volume_cyst, fill=condition)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
    ggtitle('cyst volume') + theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90))
  
  ggplot(cyst_all, aes(x = condition, y=sphericity_cyst, fill=condition)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
    ggtitle('cyst volume') + theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90))
  
  p2 = ggplot(cyst_all, aes(x = volume.log10)) +
    geom_histogram(binwidth = 0.1)
  
  p3 = ggplot(cyst_all, aes(x = sphericity_cyst)) +
    geom_histogram(binwidth = 0.01)
  
  p23 = ggplot(cyst_all, aes(x = volume.log10, y = sphericity_cyst, color = condition)) +
    geom_point(size = 1) +
    geom_hline(yintercept=0.8, colour = "red") + geom_vline(xintercept = 3.5, colour = "red")
  
  sels = which(cyst_all$volume.log10 >=3.5 & cyst_all$sphericity_cyst >=0.8)
  cyst_sel = cyst_all[sels, ]
  
  p4 = as_tibble(cyst_sel) %>% 
    group_by(condition) %>% tally() %>%
    ggplot(aes(x = condition, y = n, fill = condition)) +
    geom_bar(stat = "identity") +
    theme(legend.position = "none")  + 
    ggtitle('nb of cysts ') +
    theme(axis.text.x = element_text(angle = 90))
  
  p5 = ggplot(cyst_sel, aes(x = condition, y=AreaShape_Volume_cyst, fill=condition)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
    ggtitle('cyst volume') + theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90))
  
  ggplot(cyst_sel, aes(x = condition, y=sphericity_cyst, fill=condition)) + 
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
  
}


##########################################
# fp filtering
##########################################
res = readRDS(file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_cystFilering_', analysis.verison, '.rds'))

cond.id = paste0(res$condition, '_', res$ID_fp)

fp_all = res[match(unique(cond.id), cond.id), ]
fp_all = fp_all[!is.na(fp_all$ID_fp), ]

fp_all$volume.log10 = log10(fp_all$AreaShape_Volume_fp)
fp_all$volume.ratio = (fp_all$AreaShape_Volume_fp/fp_all$AreaShape_Volume_cyst)
fp_all$volume.ratio.log10 = log10(fp_all$AreaShape_Volume_fp/fp_all$AreaShape_Volume_cyst)
fp_all$dist.cyst.radius.ratio = fp_all$Distance_Centroid_organoid_fp/(fp_all$AreaShape_EquivalentDiameter_cyst/2)

ggplot(fp_all, aes(x = volume.log10, y = volume.ratio.log10)) +
  geom_point(size = 0.2) + 
  geom_hline(yintercept = log10(0.002), colour = "red") + 
  geom_vline(xintercept = 2, colour = "red") + 
  geom_vline(xintercept = 2.5, colour = "red") +
  ggtitle('fp/cyst volume ratio and volume in log10') 

ggplot(fp_all, aes(x = volume.log10, y=Intensity_MeanIntensity_FOXA2_fp, fill=condition)) + 
  geom_point(size = 0.2) + 
  geom_hline(yintercept = 0.01, colour = "red") + 
  geom_hline(yintercept = 0.15, colour = "red") + 
  geom_vline(xintercept = 2, colour = "red") + 
  geom_vline(xintercept = 2.5, colour = "red") + 
  ggtitle('fp mean foxa2 intensity and volume in log10') 

cat('Please specify the threshods of fp filtering \n')

sels = which(fp_all$volume.log10 > 2.5 & fp_all$Intensity_MeanIntensity_FOXA2_fp >0.01 &
               fp_all$Intensity_MeanIntensity_FOXA2_fp < 0.15 & 
               fp_all$volume.ratio <= 1.0 )

fp_sel = fp_all[sels, ]

res = res[which(is.na(res$ID_fp) | !is.na(match(res$ID_fp, fp_sel$ID_fp))), ]

saveRDS(res, file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_cyst.fp.Filering_', analysis.verison, '.rds'))


cat('Do you wannt to compare the parameters before and after cyst filtering')
Compare.fp.fitering.before.vs.after = TRUE

if(Compare.fp.fitering.before.vs.after){
  p1 = ggplot(fp_all, aes(x = condition, y=AreaShape_Volume_fp, fill=condition)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
    ggtitle('fp volume') + theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, size = 10)) +
    scale_y_continuous(trans='log10')
  
  p2 = ggplot(fp_all, aes(x = condition, y=sphericity_fp, fill=condition)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
    ggtitle('fp sphericity') + theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, size = 10))
  
  p3 = ggplot(fp_all, aes(x = condition, y=Intensity_MeanIntensity_FOXA2_fp, fill=condition)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.1) + 
    ggtitle('fp mean FoxA2 intensity') + theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, size = 10))
  
  p4 = ggplot(fp_all, aes(x = volume.log10)) +
    geom_histogram(binwidth = 0.1)
  
  p5 = ggplot(fp_all, aes(x = sphericity_fp)) +
    geom_histogram(binwidth = 0.02)
  
  p6 = ggplot(fp_all, aes(x = Intensity_MeanIntensity_FOXA2_fp)) +
    geom_histogram(binwidth = 0.001) +
    geom_vline(xintercept = 0.01, colour = "red")
  
  p45 = ggplot(fp_all, aes(x = volume.log10, y = sphericity_fp)) +
    geom_point(size = 0.2) +
    #geom_hline(yintercepst=0.85, colour = "red") + 
    geom_vline(xintercept = 1.5, colour = "red") + 
    geom_vline(xintercept = 2, colour = "red")
  
  p56 = ggplot(fp_all, aes(x = volume.log10, y = Intensity_MeanIntensity_FOXA2_fp)) +
    geom_point(size = 0.2) +
    geom_hline(yintercept=0.01, colour = "red") + 
    geom_hline(yintercept=0.15, colour = "red") + 
    geom_vline(xintercept = 1.5, colour = "red") + 
    geom_vline(xintercept = 2, colour = "red")
  
  
  p7 = ggplot(fp_sel, aes(x = condition, y=volume.log10, fill=condition)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
    ggtitle('fp volume') + theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, size = 10))
  
  p8 = ggplot(fp_sel, aes(x = condition, y=sphericity_fp, fill=condition)) + 
    geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, size = 0.5) + 
    ggtitle('fp sphericity') + theme(legend.position = "none") +
    theme(axis.text.x = element_text(angle = 90, size = 10))
  
  p9 = ggplot(fp_sel, aes(x = condition, y=Intensity_MeanIntensity_FOXA2_fp, fill=condition)) + 
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
}


########################################################
########################################################
# Section : extract and visualize the turing-model relevant parameters
#  
########################################################
########################################################
##########################################
# extract turing-relevant parameters
######################################### 
source('Functions_3Dimage.R')

res = readRDS(file = paste0(Rdata, '/mergedTable_cyst.fp_allConditions_cyst.fp.Filering_', analysis.verison, '.rds'))
cat(length(unique(res$ID_cyst)), ' cysts and ', length(unique(res$ID_fp)) -1, 'fps\n')

params = extract.turing.parameters.cellProfiler(res, pixel.scale = 3, cyst.overlapRatio.threshold = 0.01)


##########################################
# plot parameters with given conditions
##########################################
conds = unique(params$condition)
cat('all conditions found are: \n', paste0(as.character(conds), collapse = '\n'), '\n')

cat('please specify conditions you want to show \n')
#conditions_select = c("N2B27", "RA", 
#                'RA+BMP0.5', 'RA+BMP1', 'RA+BMP1.5', 'RA+BMP2', 'RA+BMP2.5', 'RA+BMP5', 'RA+BMP7.5', 'RA+BMP10')

conditions_select = c("N2B27", "RA", 
                'RA+SAG1nM', 'RA+SAG5nM', 'RA+SAG10nM', 'RA+SAG25nM', 'RA+SAG50nM', 'RA+SAG100nM', 'RA+SAG500nM', 'RA+SAG1000nM',
                'RA+Cyclo1uM', 'RA+Cyclo5uM')

pdfname = paste0(resDir, '/NTorganoid_mouse_Teresa_quickAnalysis_titration_SAG.Cyclo.pdf')

sels = which(!is.na(match(params$condition, conditions_select)))
nb.fp = as.numeric(as.character(params$nb.fp[sels]))
print(table(nb.fp))

sels = sels[which(nb.fp>=0 & nb.fp<7)]

p0 = as_tibble(params[sels, ]) %>% 
  group_by(condition, nb.fp) %>% tally() %>%
  ggplot(aes(x = factor(condition, levels = conditions_select), y = n, fill = nb.fp)) +
  geom_bar(stat = "identity") +
  theme_classic() + ggtitle('nb of cysts and fp nb distribution') +
  theme(axis.text.x = element_text(angle = 90, size = 16))


p01 = as_tibble(params[sels, ]) %>% 
  group_by(condition, nb.fp) %>% tally() %>%  mutate(percent = n/sum(n)) %>%
  ggplot(aes(x = factor(condition, levels = conditions_select), y = percent, fill = nb.fp)) +
  geom_bar(stat = "identity") +
  theme_classic() + ggtitle(' fractions of cyst with x fp ') +
  theme(axis.text.x = element_text(angle = 90, size = 16))

p1 = ggplot(params[sels, ], aes(x = factor(condition, levels = conditions_select), y=volume, fill=condition)) + 
  geom_violin(width = 1.2) + ggtitle('cyst volume') +
  theme(axis.text.x = element_text(angle = 90, size = 16))

p2 = ggplot(params[sels, ], aes(x = factor(condition, levels =conditions_select), y=overlap.ratio, fill=condition)) + 
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.1) +  ylim(0, 0.75) + 
  ggtitle('cyst fraction overlapped by fp') +
  theme(axis.text.x = element_text(angle = 90, size = 16))

p3 = ggplot(params[sels[which(as.numeric(as.character(params$nb.fp[sels]))>=1)], ], aes(x=factor(condition, levels = conditions_select), y=volume.fp, fill=condition)) + 
  geom_violin(width = 1.5) + ggtitle('foxa2 volume') +
  theme(axis.text.x = element_text(angle = 90, size = 16)) 

p4 = ggplot(params[sels[which(as.numeric(as.character(params$nb.fp[sels]))>=1)], ], aes(x = factor(condition, levels = conditions_select), y=foxa2.fp, fill=condition)) + 
  geom_violin(width = 1.5) + ggtitle('FoxA2 mean intensity') + ylim(0, 20000) + 
  theme(axis.text.x = element_text(angle = 90, size = 16))

p5 = ggplot(params[sels, ], aes(fill=condition, y=olig2 , x = condition)) + 
  geom_violin() + ggtitle('Olig2 mean intensity') +
  theme(axis.text.x = element_text(angle = 90, size = 16))

p6 = ggplot(params[sels, ], aes(x=nb.fp, y=volume, color=condition, fill = condition)) +
  geom_boxplot() + ggtitle('size dependency of fp nb (cyst volume)') + 
  theme(axis.text.x = element_text(angle = 0, size = 16), 
        axis.text.y = element_text(size = 16, angle = 90) )

p61 = ggplot(params[sels, ], aes(x=nb.fp, y=radius.cyst, color=condition, fill = condition)) +
  geom_violin() + ggtitle('size dependency of fp nb (cyst radius)') +
  theme(axis.text.x = element_text(angle = 90, size = 16))

p7 = ggplot(params[sels[which(as.numeric(as.character(params$nb.fp[sels]))>1)], ], 
            aes(x=volume, y=dist.fp, color=condition)) +
  geom_point(size = 2.5) + ggtitle('distance between fps (wavelength)') +
  theme(axis.text.x = element_text(angle = 90, size = 16))

p8 = ggplot(params[sels[which(as.numeric(as.character(params$nb.fp[sels]))>1)], ], 
            aes(x=factor(condition, levels = conditions_select), y=dist.fp, color=condition, fill = condition)) +
  #geom_violin() + 
  geom_boxplot() +
  #geom_jitter(width = 0.1, color = 'black', size = 1.0) +
  ggtitle('distance between fps (wavelength)') +
  theme(axis.text.x = element_text(angle = 90, size = 16))

#p9 = ggplot(params[which(params$condition == 'R1.normRA.d6'|params$condition == 'R1.advRA.d6' ), ], aes(x=nb.fp, y=volume, color=condition, fill = condition)) +
#  geom_boxplot() + ggtitle('size dependency of fp nb (cyst volume)') 

pdf(pdfname,  width = 18, height = 10)
grid.arrange(p0, nrow = 1, ncol = 1)
grid.arrange(p01, nrow = 1, ncol = 1)
grid.arrange(p1, nrow = 1, ncol = 1)
grid.arrange(p2, nrow = 1, ncol = 1)
grid.arrange(p5, nrow = 1, ncol = 1)
grid.arrange(p7, nrow = 1, ncol = 1)
grid.arrange(p8, nrow = 1, ncol = 1)
grid.arrange(p4, nrow = 1, ncol = 1)
grid.arrange(p3, nrow = 1, ncol = 1)
grid.arrange(p6, nrow = 1, ncol = 1)
grid.arrange(p61, nrow = 1, ncol = 1)
#grid.arrange(p9, nrow = 1, ncol = 1)

dev.off()

##########################################
# in case need to save result in table
##########################################

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

