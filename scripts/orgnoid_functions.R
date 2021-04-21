##########################################################################
##########################################################################
# Project: Orgnoid patterning 
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Apr 16 15:15:20 2021
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
    kk = which(param[, 1] != ' ')
    
    name.param = gsub('.csv', '', name.param)
    name.param = gsub(' ', '.', name.param)
    name.param = gsub('=', '', name.param)
    name.param = gsub('Surfaces_[1-9]_', '', name.param)
    name.param = gsub('_', '.', name.param)
    
    cat(conds[n],  ' : m = ', m, ' param.name :',  name.param, '\n')
    
    # colnames of extracted parameterse
    # I hate coding like 
    names = as.character(unlist(param[kk[2], ]))
    names = gsub('_', '.', names)
    names = gsub(' ', '.', names)
    
    param = data.frame(param[(kk[2]+1):nrow(param), ])
    colnames(param) = names
    index.id = which(names == 'ID')
    
    if(length(index.id) != 1) {
      stop('Error : no ID column found !!! \n')
      
    }else{
      # keep the original image name and id
      if(m == 1){
        images.names = param[, c(which(colnames(param) == 'Original.Image.Name'), 
                                 which(colnames(param) == 'Original.Image.ID'))]
      }
      
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
          
        }#else{
        ## possible a bug here
        #param = param[, c(1:index.unit)]
        #}
      }
      
      colnames(param)[2:ncol(param)] = paste0(colnames(param)[2:ncol(param)], '_', name.param)
      
      if(m == 1 ){
        #colnames(param)[2:ncol(param)] = paste0(colnames(param)[2:ncol(param)], '_', name.param)
        res.cc = data.frame(images.names, param, stringsAsFactors = FALSE)
        
        
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

verify.overlapping.by.distance = function(x0, yy)
{
  # x0 = res.cyst[n,]; yy = res.fp[kk,]
  rads = (as.numeric(x0$Volume.Unit._Volume)/(4/3*pi))^(1/3)
  dists = c()
  for(m in 1:nrow(yy))
  {
    y0 = yy[m,]
    dists = c(dists, sqrt((as.numeric(x0$Position.X.Unit.µm_Position) - as.numeric(y0$Position.X.Unit.µm_Position))^2 +
                          (as.numeric(x0$Position.Y.Unit.µm_Position) - as.numeric(y0$Position.Y.Unit.µm_Position))^2 +
                          (as.numeric(x0$Position.Z.Unit.µm_Position) - as.numeric(y0$Position.Z.Unit.µm_Position))^2))
  }
  
  overlapped = dists <= rads
  return(overlapped)
  
}


##########################################
# function that serach all fp for each cyst
# but this logic doesn't work quite well, because we don't know exactly the radius of cyst due to the missing lumen volumes in the 
# segmentation
##########################################
find.fp.for.each.cyst = function(res.cyst, res.fp)
{
  # res.cyst = res1; res.fp = res2
  res.cyst = data.frame(res.cyst, stringsAsFactors = FALSE)
  res.fp = data.frame(res.fp, stringsAsFactors = FALSE)
  
  index.cyst = c()
  index.fp = c()
  
  for(n in 1:nrow(res.cyst)){
    # n = 11
    volume.ratio.overlapped = as.numeric(res.cyst[n, grep('Overlapped.Volume.Ratio.to.Surfaces', colnames(res.cyst))])
    shortest.dist = as.numeric(res.cyst[n, grep('Shortest.Distance.to.Surfaces', colnames(res.cyst))])
    
    # overlapped by fp indicated by overlapped ratio and shortest distance between surfaces
    if(volume.ratio.overlapped > 0.01 & shortest.dist < 10^-2){
      
      # search the fp in the same image 
      kk = which(res.fp$Original.Image.Name == res.cyst$Original.Image.Name[n])
      volume.ratio.overlapped.fp = as.numeric(res.fp[kk, grep('Overlapped.Volume.Ratio.to.Surfaces', colnames(res.fp))])
      shortest.dist.fp = as.numeric(res.fp[kk, grep('Shortest.Distance.to.Surfaces', colnames(res.fp))])
      
      kk = kk[which(volume.ratio.overlapped.fp > 0.95 & shortest.dist.fp < 10^-2)]
      
      op = verify.overlapping.by.distance(res.cyst[n,], res.fp[kk, ])
      kk = unique(kk[op])
      
      if(length(kk) ==0){
        cat('0 FP in this cyst -- \n')
        index.cyst = c(index.cyst, n)
        index.fp = c(index.fp, NA)
      }else{
        cat(length(kk),  ' FP found \n')
        index.cyst = c(index.cyst, rep(n, length(kk)))
        index.fp = c(index.fp, kk)
      }
       
    }else{ # no overlapping volum and shortest distance is larger than 0
     cat('0 FP in this cyst --  \n')
      index.cyst = c(index.cyst, n)
      index.fp = c(index.fp, NA)
    }
    
  }
  
  colnames(res.cyst) = paste0(colnames(res.cyst), '_cyst')
  colnames(res.fp) = paste0(colnames(res.fp), '_fp')
  
  res.cp = data.frame(res.cyst[index.cyst, ], res.fp[index.fp, ], stringsAsFactors = FALSE)
  
  return(res.cp)
  
}

##########################################
# assign fp to cyst by calculate the distances between centers
##########################################
calculate.distance = function(x0, yy)
{
  # x0 = res.cyst[n,]; yy = res.fp[kk,]
  dists = c()
  for(m in 1:nrow(yy))
  {
    y0 = yy[m,]
    dists = c(dists, sqrt((as.numeric(x0$Position.X.Unit.µm_Position) - as.numeric(y0$Position.X.Unit.µm_Position))^2 +
                            (as.numeric(x0$Position.Y.Unit.µm_Position) - as.numeric(y0$Position.Y.Unit.µm_Position))^2 +
                            (as.numeric(x0$Position.Z.Unit.µm_Position) - as.numeric(y0$Position.Z.Unit.µm_Position))^2))
  }
  
  return(dists)
  
}

find.closest.cyst = function(x0, yy)
{
  # x0 = res.fp[n,]; yy = res.cyst[kk,]
  #rads = (as.numeric(x0$Volume.Unit._Volume)/(4/3*pi))^(1/3)
  dists = c()
  for(m in 1:nrow(yy))
  {
    y0 = yy[m, ]
    dists = c(dists, sqrt((as.numeric(x0$Position.X.Unit.µm_Position) - as.numeric(y0$Position.X.Unit.µm_Position))^2 +
                            (as.numeric(x0$Position.Y.Unit.µm_Position) - as.numeric(y0$Position.Y.Unit.µm_Position))^2 +
                            (as.numeric(x0$Position.Z.Unit.µm_Position) - as.numeric(y0$Position.Z.Unit.µm_Position))^2))
  }
  
  return(which(dists == min(dists))[1])
  
}

fine.correction.with.dist.volume = function(mapping, res.cyst, res.fp)
{
  for(n in 1:length(mapping))
  {
    # n = 98
    if(length(mapping[[n]])>0){
      index.cyst = which(res.cyst$ID == gsub('cyst_','', names(mapping)[n]))
      index.fp = mapping[[n]]
      
      # verfiy volume
      overlap.volumes.cyst = as.numeric(res.cyst[index.cyst, grep('Overlapped.Volume.to.Surfaces', colnames(res.cyst))])
      overlap.volume.fps = as.numeric(res.fp[index.fp, grep('Overlapped.Volume.to.Surfaces', colnames(res.fp))])
      volume.fps = as.numeric(res.cyst[index.fp, grep('Volume.Unit', colnames(res.fp))])
      
      ## verify the distance
      dists = calculate.distance(res.cyst[index.cyst, ], res.fp[index.fp, ])
      if(max(dists)/min(dists) > 5){
        cat(n, '--', names(mapping)[n], 'to check -- dist: ')
        print(dists)
      }
      
      # test the furthest fp
      Test.new.cyst = FALSE
      if(Test.new.cyst){
        ii = index.fp[which(dists == max(dists))]
        jj = which(res.cyst$Original.Image.Name == res.fp$Original.Image.Name[ii])
        ds = calculate.distance(res.fp[ii, ], res.cyst[jj, ])
        
        jj = jj[order(ds)]
        ds = ds[order(ds)]    
      }
    }
  }
  
}

find.cyst.for.each.fp = function(res.cyst, res.fp)
{
  # res.cyst = res1; res.fp = res2
  res.cyst = data.frame(res.cyst, stringsAsFactors = FALSE)
  res.fp = data.frame(res.fp, stringsAsFactors = FALSE)
  
  res.cyst$ID = gsub(' ', '', res.cyst$ID)
  res.fp$ID = gsub(' ', '', res.fp$ID)
  
  # define a list in which each element named with cyst_ID and its contents will be a vector of assocaited fp IDs if there are any
  mapping = vector("list", length = length(unique(res.cyst$ID)))
  names(mapping) = paste0('cyst_', unique(res.cyst$ID))
  
  # loop over fp to serach for associated cyst
  for(n in 1:nrow(res.fp)){
    # n = 1
    cat('index of fp : ', n, '\n')
    vratio.overlapped = as.numeric(res.fp[n, grep('Overlapped.Volume.Ratio.to.Surfaces', colnames(res.fp))])
    #shortest.dist = as.numeric(res.fp[n, grep('Shortest.Distance.to.Surfaces', colnames(res.fp))])
    
    # proceed if overlapped volume > 0.01 minimum
    if(vratio.overlapped > 0.01){
      # search cysts in the same image as fp and select those with overlapped volumes 
      kk = which(res.cyst$Original.Image.Name == res.fp$Original.Image.Name[n])
      kk = kk[which(as.numeric(res.cyst[kk, grep('Overlapped.Volume.Ratio.to.Surfaces', colnames(res.cyst))]) > 0)]
      
      if(length(kk) > 1){
        kk.closest = kk[find.closest.cyst(res.fp[n,], res.cyst[kk, ])]
        mapping[[kk.closest]] = c(mapping[[kk.closest]], n)
      }else{
        if(length(kk) == 1){
          kk.closest = kk
          mapping[[kk.closest]] = c(mapping[[kk.closest]], n)
        } 
      }
    }
  }
  
  #mapping.c = fine.correction.with.dist.volume(mapping, res.cyst, res.fp)
  
  nb.fp = c()
  index.cf = c()
  for(n in 1:length(mapping))
  {
    nb.fp = c(nb.fp, length(mapping[[n]]))
    index.cyst = which(res.cyst$ID == gsub('cyst_','', names(mapping)[n]))
    index.fp = mapping[[n]]
    
    index.cf = rbind(index.cf, cbind(rep(index.cyst, length(index.fp)), index.fp))
  }
  hist(nb.fp, breaks = c(-1:max(nb.fp)))
  
  colnames(res.cyst) = paste0(colnames(res.cyst), '_cyst')
  colnames(res.fp) = paste0(colnames(res.fp), '_fp')
  
  res.cp = data.frame(res.cyst[index.cf[,1], ], res.fp[index.cf[,2], ], stringsAsFactors = FALSE)
  
  return(res.cp)
  
  
}


##########################################
# extract turing model related parameters 
##########################################
calcuate.cyst.radius = function(xx0)
{
  # xx0 = res.cp[kk, ]
  
  # first check if cyst id is unique
  if(length(unique(xx0$ID_cyst) == 1)){
    
    x0 = xx0[, grep('_cyst$', colnames(xx0))]
    yy = xx0[, grep('_fp$', colnames(xx0))]
    colnames(x0) = gsub('_cyst$', '', colnames(x0))
    colnames(yy) = gsub('_fp$', '', colnames(yy))
    dd = calculate.distance(x0[1, ], yy)
    rfp = (as.numeric(yy$Volume.Unit._Volume)/(4/3*pi))^(1/3)
    rrc = median(dd + rfp)
    return(rrc)
    
  }else{
    stop('multiple cyst ID found !!!')
  }
}

extract.turing.parameters = function(res.cp)
{
  # res.cp = res
  nb.fp = c()
  radius = c()
  
  cyst.id = unique(res.cp$ID_cyst)
  
  for(n in 1:length(cyst.id))
  {
    kk = which(res.cp$ID_cyst == cyst.id[n])
    if(length(kk) == 0){
      nb.fp = c(nb.fp, 0)
      radius = c(radius, NA)
    }
    
    if(length(kk) >=1){
      nb.fp = c(nb.fp, length(kk))
      radius = c(radius, calcuate.cyst.radius(res.cp[kk, ]))
    }
        
  }
  
  
}


