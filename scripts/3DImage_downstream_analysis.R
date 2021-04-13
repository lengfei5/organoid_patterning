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

# specific input and output folders
dataDir = '../outputs/210217_hNTdrugs3_0310outputs'
resDir = 'results'
analysis.verison = '_test20210410'

cyst.channel = '1'
floorplat.channel = '2'

if(!dir.exists(resDir)) dir.create(resDir)
########################################################
########################################################
# Section : concatinate the readouts of 3D image analysis
#  
########################################################
########################################################
conditions.list = dir(path = dataDir, pattern = '*_Statistics', full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
conds = basename(conditions.list)
conds = gsub('_Statistics', '', conds)

res = c()

for(n in 1:length(conds))
{
  # n = 1
  cat(n, ' -- start ', conds[n], '\n')
  files.cyst = list.files(path = paste0(conditions.list[n], '/Surfaces_', cyst.channel,   '_Statistics'), 
                          pattern = '*.csv', full.names = TRUE)
  
  
}



