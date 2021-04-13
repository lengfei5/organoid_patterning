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

rm(list=ls())

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

##########################################
# functions 
##########################################
cat.image.parameters.per.condition = function(file.surface)
{
  name.param = basename(files.cyst[m])
  #cat(name.param, ' -- ')
  param = read.csv(file = files.cyst[m], sep = ',', header = FALSE, comment.char = '=', stringsAsFactors = FALSE)
  kk = which(param[,1] != ' ')
  
  name.param = gsub('.csv', '', name.param)
  name.param = gsub(' ', '.', name.param)
  name.param = gsub('=', '', name.param)
  cat(conds[n],  'm = ', m, ' param.name :',  name.param, '\n')
  
  names = as.character(unlist(param[kk[2], ]))
  param = data.frame(param[(kk[2]+1):nrow(param), ])
  colnames(param) = names
  index.id = which(names == 'ID')
  
  if(length(index.id) != 1) {
    cat('Error : no ID column found !!! \n')
  }else{
    param = param[, c(index.id, 1:(index.id-1))]
    index.unit = grep('Unit', colnames(param))
    if(length(index.unit) >= 1){
      index.unit = max(index.unit)
      param = param[, c(1:index.unit)]
    }
    
    colnames(param)[1:ncol(param)] = paste0(colnames(param)[1:ncol(param)], '_', name.param)
    
  }
}


##########################################
# 
##########################################
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
  
  res.cc = c()
  for(m in 1:length(files.cyst))
  {
    # m = 1
    if(m == 1){
      colnames(param)[2:ncol(param)] = paste0(colnames(param)[2:ncol(param)], '_', name.param)
      res.cc = data.frame(param, stringsAsFactors = FALSE)
      
    }else{
      
      ii = match(res.cc$ID, param$ID)
      res.cc = data.frame(res.cc, param[ii, ])  
    }
    
  }
  
}



