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
analysis.verison = 'hNTdrugs3_0310_Analysis20210410'

cyst.channel = '1'
floorplat.channel = '2'

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(tabDir)) dir.create(tabDir)

source('orgnoid_functions.R')

save.table.each.condition = FALSE

##########################################
# find associated fp for cyst at each condition
##########################################
conditions.list = dir(path = dataDir, pattern = '*_Statistics', full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
conds = basename(conditions.list)
conds = gsub('_Statistics', '', conds)

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
  
  if(save.table.each.condition){
    write.table(res, file = paste0(tabDir, '/', analysis.verison, '_condition_', conds[n], '_cyst_fp.txt'), 
                sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
  }
  
}

########################################################
########################################################
# Section : extract relevant parameters and compare across conditions
# 
########################################################
########################################################


# filter cyst or/and floorplates using global parameters

# extract turing-relevant parameters
source('orgnoid_functions.R')

params = extract.turing.parameters(res)





