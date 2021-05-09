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
dataDir = '../210217_hNTdrugs3_0310outputs'
resDir = '../results'
tabDir = paste0(resDir, '/tables')
Rdata = paste0(resDir, '/Rdata')
analysis.verison = 'hNTdrugs3_0310_Analysis20210507'

cyst.channel = '1'
floorplat.channel = '2'

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(tabDir)) dir.create(tabDir)
if(!dir.exists(dataDir)) dir.create(dataDir)

save.table.each.condition = FALSE

##########################################
# find associated fp for cyst at each condition
##########################################
conditions.list = dir(path = dataDir, pattern = '*_Statistics', full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
conds = basename(conditions.list)
conds = gsub('_Statistics', '', conds)

source('orgnoid_functions.R')

for(n in 1:length(conds))
{
  # n = 18
  cat(n, ' -- start ', conds[n], '\n')
  
  files.cyst = list.files(path = paste0(conditions.list[n], '/Surfaces_', cyst.channel,   '_Statistics'), 
                          pattern = '*.csv', full.names = TRUE)
  files.cyst =  files.cyst[grep('_Overall', files.cyst, invert = TRUE)]
  
  res1 = cat.image.parameters.per.condition(files.cyst)
  
  files.fp = list.files(path = paste0(conditions.list[n], '/Surfaces_', floorplat.channel,   '_Statistics'), 
                        pattern = '*.csv', full.names = TRUE)
  files.fp =  files.fp[grep('_Overall', files.fp, invert = TRUE)]
  
  res2 = cat.image.parameters.per.condition(files.fp)
  
  if(save.table.each.condition){
    write.table(res1, file = paste0(resDir, '/', analysis.verison, '_condition_', conds[n], '_cyst.txt'), 
                sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
    write.table(res2, file = paste0(resDir, '/', analysis.verison, '_condition_', conds[n], '_floorplate.txt'), 
                sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
  }
  
  
  source('orgnoid_functions.R')
  
  res = find.cyst.for.each.fp(res.cyst = res1, res.fp = res2)
  
  write.table(res, file = paste0(tabDir, '/', analysis.verison, '_condition_', conds[n], '_cyst_fp.txt'), 
                sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  
}


# merge tables of all conditions
files = list.files(path = tabDir, pattern = '_cyst_fp.txt', full.names = TRUE)
files = files[grep(analysis.verison, files)]

cc = c()
res = c()
for(n in 1:length(files))
{
  # n = 1
  cat(n, ' \t')
  xx = read.table(files[n], sep = '\t', header = TRUE)
  colnames(xx)[grep('Overlapped.Volume.Ratio', colnames(xx))] = c('Overlapped.Volume.Ratio_cyst', 'Overlapped.Volume.Ratio_fp')
  colnames(xx)[grep('Overlapped.Volume.to.Surfaces.Unit', colnames(xx))] = c('Overlapped.Volume.to.Surfaces.Unit_cyst', 
                                                                             'Overlapped.Volume.to.Surfaces.Unit_fp')
  colnames(xx)[grep('Shortest.Distance.to.Surfaces', colnames(xx))] = c('Shortest.Distance.to.Surfaces_cyst', 
                                                                        'Shortest.Distance.to.Surfaces_fp')
  
  c = gsub(paste0(analysis.verison, '_condition_'),'', basename(files[n]))
  c = gsub('_cyst_fp.txt', '', c)
  
  if(n == 1){
    res = xx
  }else{
    res = rbind(res, xx)
    
  }
  cc = c(cc, rep(c, nrow(xx)))
  cat(nrow(xx),  'row,  done \n')
}

res = data.frame(condition = cc, res, stringsAsFactors = FALSE)

saveRDS(res, file = paste0(Rdata, 'mergedTable_cyst.fp_allConditions_', analysis.verison, '.rds'))

Double.check = FALSE
if(Double.check){
  res = readRDS(file = paste0(Rdata, 'mergedTable_cyst.fp_allConditions.rds'))
  
  jj = which(!is.na(res$ID_fp))
  
  head(res[jj, grep('Original.Image.Name', colnames(res))])
  
  length(which(as.character(res$Original.Image.Name_cyst[jj]) != as.character(res$Original.Image.Name_fp[jj])))
  
  res = res[which(res$condition == 'RA_LDNSB' & res$Original.Image.Name_cyst == '210217_nodrug_LDNSB_4_01_[ims1_2021-03-09T16-53-43.863]'),
            ]
  
  xx = (cbind(res$OriginalID_cyst, res$OriginalID_fp, res$Distance.to.Image.Border.XY.Unit.µm_Distance.to.Image.Border.XY.Img1_cyst))
  colnames(xx) = c('cyst.ID', 'fp.id', 'Cyst.distanct.image.board.XY')
  
  write.table(xx, file = paste0(resDir, '/table_for_Hannah_manualCheck_cyst_fp.assignment_refined.Assignement.txt'),
              sep = '\t', col.names = TRUE, row.names = FALSE, quote = FALSE)
  
  
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


