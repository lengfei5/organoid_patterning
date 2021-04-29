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
# Section : concatinate the readouts of 3D image analysis
# map the floorplate with the cysts 
########################################################
########################################################
rm(list=ls())

# specific input and output folders
dataDir = '../210217_hNTdrugs3_0310outputs'
resDir = '../results'
tabDir = paste0(resDir, '/tables')
Rdata = paste0(resDir, '/Rdata')
analysis.verison = 'hNTdrugs3_0310_Analysis20210410'

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
  # n = 12
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

saveRDS(res, file = paste0(Rdata, 'mergedTable_cyst.fp_allConditions.rds'))

########################################################
########################################################
# Section : extract relevant parameters and compare across conditions
# 
########################################################
########################################################
res = readRDS(file = paste0(Rdata, 'mergedTable_cyst.fp_allConditions.rds'))

##########################################
# filter cyst or/and floorplates using global parameters
##########################################
res = res[which(res$Distance.to.Image.Border.XY.Unit.µm_Distance.to.Image.Border.XY.Img1_cyst > 0), ]

res = res[which(res$Overlapped.Volume.Ratio_fp > 0.25 | is.na(res$Overlapped.Volume.Ratio_fp)), ]

res = res[which(res$Sphericity.Unit._Sphericity_cyst > 0.8), ]

##########################################
# extract turing-relevant parameters
#########################################

source('orgnoid_functions.R')

params = extract.turing.parameters(res)



##########################################
# visualize the turing-model relevant parameters
##########################################
require(ggplot2)
require(grid)
require(gridExtra)
params$condition = as.factor(params$condition)
params$nb.fp = as.factor(params$nb.fp)
params$radius.cyst = as.numeric(params$radius.cyst)
params$volume = as.numeric(params$volume)
params$dist.fp = as.numeric(params$dist.fp)
params$foxa2.fp = as.numeric(params$foxa2.fp)
params$radius.fp = as.numeric(params$radius.fp)

conds = unique(params$condition)


# check controls first
conds.sels = list(
  # controls 
  which(params$condition == "RA_LDNSB"|
              params$condition == 'noRA_LDNSB'| 
              params$condition == 'RA_noLDNnoSB'|
              params$condition == 'RA_LDNonly'),
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
          params$condition == 'RA_SBonly') 

)

pdfname = paste0(resDir, '/perturbation_summary.pdf')
pdf(pdfname,  width = 24, height = 16)

for(n in 1:length(conds.sels))
{
  # n = 1
  sels = conds.sels[[n]]
  nb.fp = as.numeric(as.character(params$nb.fp[sels]))
  sels = sels[which(nb.fp>=0 & nb.fp<10)]
  
  p1 = ggplot(params[sels, ], aes(x=nb.fp, y=radius.cyst, color=condition, fill = condition)) +
    geom_violin() + ggtitle('estimated cyst radius')
  
  
  p2 = ggplot(params[sels, ], aes(x=nb.fp, y=dist.fp, fill=condition)) +
    geom_violin() + ggtitle('distance between fps')
  
  
  p3 = ggplot(params[sels, ], aes(x=nb.fp, y=foxa2.fp, fill=condition)) +
    geom_violin() + ggtitle('foxa2 mean intensity')
  
  p4 = ggplot(params[sels, ], aes(x=nb.fp, y=radius.fp, fill=condition)) +
    geom_violin() + ggtitle('foxa2 radius')
  
  grid.arrange(p1, p2, p3, p4, nrow = 2, ncol = 2)
  
}

dev.off()




