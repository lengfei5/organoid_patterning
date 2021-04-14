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

##########################################
# functions 
##########################################
cat.image.parameters.per.condition = function(files.surface)
{
  # files.surface = files.cyst;
  res.cc = c()
    
  for(m in 1:length(files.surface))
  {
    # m = 1
    name.param = basename(files.surface[m])
    #cat(name.param, ' -- ')
    param = read.csv(file = files.surface[m], sep = ',', header = FALSE, comment.char = '=', stringsAsFactors = FALSE)
    kk = which(param[,1] != ' ')
    
    name.param = gsub('.csv', '', name.param)
    name.param = gsub(' ', '.', name.param)
    name.param = gsub('=', '', name.param)
    name.param = gsub('Surfaces_[1-9]_', '', name.param)
    name.param = gsub('_', '.', name.param)
    
    cat(conds[n],  ' : m = ', m, ' param.name :',  name.param, '\n')
    
    # colnames of extracted parameterse
    names = as.character(unlist(param[kk[2], ]))
    names = gsub('_', '.', names)
    names = gsub(' ', '.', names)
    
    param = data.frame(param[(kk[2]+1):nrow(param), ])
    colnames(param) = names
    index.id = which(names == 'ID')
    
    if(length(index.id) != 1) {
      stop('Error : no ID column found !!! \n')
      
    }else{
      # select columns before the ID column
      param = param[, c(index.id, 1:(index.id-1))]
      
      # select again the columns before Unit if it exits in the table
      index.unit = grep('Unit', colnames(param))
     
      if(length(index.unit) >= 1){
        index.unit = max(index.unit)
        unit.uniq = unique(as.character(param[, index.unit]))
        unit.uniq = gsub('[\\^]', '', unit.uniq)
        if(length(unit.uniq) == 1) {
          param = param[, c(1:(index.unit-1))]
          
          colnames(param)[2:ncol(param)] = paste0(colnames(param)[2:ncol(param)], '.Unit.', unit.uniq)
          
        }else{
          param = param[, c(1:index.unit)]
        }
      }
      
      colnames(param)[2:ncol(param)] = paste0(colnames(param)[2:ncol(param)], '_', name.param)
      
      if(m == 1 ){
        #colnames(param)[2:ncol(param)] = paste0(colnames(param)[2:ncol(param)], '_', name.param)
        res.cc = data.frame(param, stringsAsFactors = FALSE)
        
      }else{
        # double check the ID match
        if(length(intersect(res.cc$ID, param$ID)) != length(res.cc$ID)) { cat('ID does not match exactly \n') }
        ii = match(res.cc$ID, param$ID)
        if(ncol(param) > 2){
          res.cc = data.frame(res.cc, param[ii, -1], stringsAsFactors = FALSE)  
        }
        if(ncol(param) == 2){
          res.cc = data.frame(res.cc, param[ii, -1], stringsAsFactors = FALSE)
          colnames(res.cc)[ncol(res.cc)] = colnames(param)[2]
        }
      }
      
    }
    
  }
  
  return(res.cc)
  
}


########################################################
########################################################
# Section : concatinate the readouts of 3D image analysis
#  
########################################################
########################################################
# rm(list=ls())

# specific input and output folders
dataDir = '../210217_hNTdrugs3_0310outputs'
resDir = '../results'
analysis.verison = 'hNTdrugs3_0310_Analysis20210410'

cyst.channel = '1'
floorplat.channel = '2'

if(!dir.exists(resDir)) dir.create(resDir)


##########################################
# 
##########################################
conditions.list = dir(path = dataDir, pattern = '*_Statistics', full.names = TRUE, recursive = FALSE, include.dirs = TRUE)
conds = basename(conditions.list)
conds = gsub('_Statistics', '', conds)


for(n in 1:length(conds))
{
  # n = 1
  cat(n, ' -- start ', conds[n], '\n')
  
  files.cyst = list.files(path = paste0(conditions.list[n], '/Surfaces_', cyst.channel,   '_Statistics'), 
                          pattern = '*.csv', full.names = TRUE)
  files.cyst =  files.cyst[grep('_Overall', files.cyst, invert = TRUE)]
  
  res1 = cat.image.parameters.per.condition(files.cyst)
  
  files.fp = list.files(path = paste0(conditions.list[n], '/Surfaces_', floorplat.channel,   '_Statistics'), 
                        pattern = '*.csv', full.names = TRUE)
  files.fp =  files.fp[grep('_Overall', files.fp, invert = TRUE)]
  
  res2 = cat.image.parameters.per.condition(files.fp)
  
  write.table(res1, file = paste0(resDir, '/', analysis.verison, '_condition_', conds[n], '_cyst.txt'), 
              sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
  write.table(res2, file = paste0(resDir, '/', analysis.verison, '_condition_', conds[n], '_floorplate.txt'), 
              sep = '\t', quote = FALSE, col.names = TRUE, row.names = FALSE)
  
  
}



