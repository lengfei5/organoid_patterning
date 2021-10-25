##########################################################################
##########################################################################
# Project: Orgnoid patterning 
# Script purpose:
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Apr 16 15:15:20 2021
##########################################################################
##########################################################################

########################################################
########################################################
# Section : process the segementation results from cellProfiler 
# 
########################################################
########################################################
merge_image.cyst.fp_fromCellProfiler = function(image, cyst, fp)
{
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
    if(n %%100 == 0) cat(n, '\n')
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
  
  return(res)
}

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
                                 which(colnames(param) == 'Original.Image.ID'), 
                                 which(colnames(param) ==  "OriginalID"))]
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
# but this logic doesn't work quite well, because we don't know exactly the radius of cyst due to 
# the missing lumen volumes in the 
# segmentation (NOT Use at the end)
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
# Assign fp to cyst by calculate the distances between centers
# Used at the end
##########################################
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

## initial motivation was to improve the assginment of fp 
## now as quality control 
fp.assigment.correction.with.dist.volume = function(mapping, res.cyst, res.fp) # not used for the moment
{
  #mappinp.saved = mapping
  newmapping = mapping
  nb.iteration = 6
  index.fp.corrected = c()
  
  for(nb.i in 1:nb.iteration)
  {
    dds = c()
    rrs = c()
    cat(' iteration of correction -- ', nb.i, '\n')
    cat(' --------------------------------------\n')
    for(n in 1:length(newmapping))
    {
      # n = 128
      if(length(newmapping[[n]])>0){
        
        index.cyst = which(res.cyst$ID == gsub('cyst_','', names(newmapping)[n]))
        index.fp = newmapping[[n]]
        
        # volume: compare the sum of fp volumes with the assigned cyst volume overlapped by fps
        overlap.volumes.cyst = as.numeric(res.cyst[index.cyst, grep('Overlapped.Volume.to.Surfaces', colnames(res.cyst))])
        overlap.volume.fps = as.numeric(res.fp[index.fp, grep('Overlapped.Volume.to.Surfaces', colnames(res.fp))])
        volume.fps = as.numeric(res.fp[index.fp, grep('Volume.Unit', colnames(res.fp))])
        rr.fp.cyst = volume.fps/overlap.volumes.cyst
        
        ## verify the distance
        dists = calculate.distance(res.cyst[index.cyst, ], res.fp[index.fp, ])
        
        dds = c(dds, dists)
        rrs = c(rrs, rr.fp.cyst)
        
        # try to find the wrongly assigned fp (furthest one here) and assign it to the second close cyst 
        if(sum(rr.fp.cyst) > 4 | max(dists)/min(dists) > 5 ){
          
          ii.fp = index.fp[which(dists == max(dists))]
          cat(n, '--', names(newmapping)[n], '-- fp index', ii.fp,  'reassigned to ')
          
          if(length(which(index.fp.corrected == ii.fp)) >= 2){
            cat('alread reassigned twice, now drop this fp \n')
            newmapping[[index.cyst]] = setdiff(newmapping[[index.cyst]], ii.fp)
            
          }else{
            index.fp.corrected = c(index.fp.corrected, ii.fp)
            ## find the second closest cyst
            kk = which(res.cyst$Original.Image.Name == res.fp$Original.Image.Name[ii.fp])
            kk = kk[which(as.numeric(res.cyst[kk, grep('Overlapped.Volume.Ratio.to.Surfaces', colnames(res.cyst))]) > 0)]
            kk = setdiff(kk, index.cyst) # exclude the cyst previously assigned 
            
            # find the closest cyst for the floor plate
            if(length(kk) > 1){
              kk.closest = kk[find.closest.cyst(res.fp[ii.fp,], res.cyst[kk, ])]
              index.newmapping = which(names(newmapping) == paste0('cyst_', res.cyst$ID[kk.closest]))
              # reassign fp to the second closet cyst and remove it from the closest one
              newmapping[[index.newmapping]] = c(newmapping[[index.newmapping]], ii.fp)
              newmapping[[index.cyst]] = setdiff(newmapping[[index.cyst]], ii.fp)
              
            }else{
              if(length(kk) == 1){
                kk.closest = kk
                index.newmapping = which(names(newmapping) == paste0('cyst_', res.cyst$ID[kk.closest]))
                # reassign fp to the second closet cyst and remove it from the closest one
                newmapping[[index.newmapping]] = c(newmapping[[index.newmapping]], ii.fp)
                newmapping[[index.cyst]] = setdiff(newmapping[[index.cyst]], ii.fp)
                
              } 
            }
            cat(kk.closest, '\n')
          }
          
          #print(dists)
          #cat('ratios of fp vs. overlapped volume of cyst')
          #print (rr.fp.cyst)
          #print(res.cyst$Sphericity.Unit._Sphericity[index.cyst])
          #cat('\n')
          #catr()
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
  
  return(newmapping)
  
}

##########################################
# Main function for fp assginment to parent cysts 
# double checked with Hannah by comparing her manual analysis
# Very low error rate were found, implying that assign the fp to the closest 
# cyst after filtering not likely cysts is very close to the true one in general
#
##########################################
find.cyst.for.each.fp = function(res.cyst, res.fp, fp.assignment.correction = TRUE, Quality.test = FALSE)
{
  # res.cyst = res1; res.fp = res2
  res.cyst = data.frame(res.cyst, stringsAsFactors = FALSE)
  res.fp = data.frame(res.fp, stringsAsFactors = FALSE)
  
  # clean the name by removing the space
  res.cyst$ID = gsub(' ', '', res.cyst$ID)
  res.fp$ID = gsub(' ', '', res.fp$ID)
  
  # define a list in which each element named with cyst_ID (the same order as res.cyst) 
  # and its contents will be a vector of assocaited fp IDs if there are any
  if(length(res.cyst$ID) != length(unique(res.cyst$ID))){
    stop('cyst ID is not unqiue')
  }else{
    #res.cyst = res.cyst[match(unique(res.cyst$ID), res.cyst)]
    mapping = vector("list", length = length(res.cyst$ID)) # d
    names(mapping) = paste0('cyst_', res.cyst$ID)
  }
  
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
      
      # find the closest cyst for the floor plate
      if(length(kk) > 1){
        kk.closest = kk[find.closest.cyst(res.fp[n,], res.cyst[kk, ])]
        index.mapping = which(names(mapping) == paste0('cyst_', res.cyst$ID[kk.closest]))
        mapping[[index.mapping]] = c(mapping[[index.mapping]], n)
      }else{
        if(length(kk) == 1){
          kk.closest = kk
          index.mapping = which(names(mapping) == paste0('cyst_', res.cyst$ID[kk.closest]))
          mapping[[index.mapping]] = c(mapping[[index.mapping]], n)
        } 
      }
    }
  }
  
  # double check the quality 
  if(fp.assignment.correction){
    mapping = fp.assigment.correction.with.dist.volume(mapping, res.cyst, res.fp)
  }
  
  #mapping.c = fine.correction.with.dist.volume(mapping, res.cyst, res.fp)
  nb.fp = c()
  index.cf = c()
  for(n in 1:length(mapping))
  {
    nb.fp = c(nb.fp, length(mapping[[n]]))
    
    index.cyst = which(res.cyst$ID == gsub('cyst_','', names(mapping)[n]))
    index.fp = mapping[[n]]
    if(length(index.fp) == 0){
      index.cf = rbind(index.cf, cbind(index.cyst, NA))
    }else{
      index.cf = rbind(index.cf, cbind(rep(index.cyst, length(index.fp)), index.fp))
    }
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

calculate.distance = function(x0, yy)
{
  # x0 = res.cyst[n,]; yy = res.fp[kk,]
  dists = c()
  for(m in 1:nrow(yy))
  {
    y0 = yy[m,]
    dists = c(dists, sqrt((as.numeric(x0$AreaShape_Center_X) - as.numeric(y0$AreaShape_Center_X))^2 +
                            (as.numeric(x0$AreaShape_Center_Y) - as.numeric(y0$AreaShape_Center_Y))^2 +
                            (as.numeric(x0$AreaShape_Center_Z) - as.numeric(y0$AreaShape_Center_Z))^2))
  }
  
  return(dists)
  
}


calculate.angle.between.two.fps = function(x, y)
{
  # x = x0[1,]; y = yy
  if(nrow(y) == 2){
    C = c(x$AreaShape_Center_X, x$AreaShape_Center_Y, x$AreaShape_Center_Y)
    A = c(y$AreaShape_Center_X[1], y$AreaShape_Center_Y[1], y$AreaShape_Center_Y[1])
    B = c(y$AreaShape_Center_X[2], y$AreaShape_Center_Y[2], y$AreaShape_Center_Y[2])
    
    a2 = sum((C-A)^2)
    b2 = sum((C-B)^2)
    c2 = sum((A-B)^2)
    alpha = acos((a2+b2 -c2)/(2*sqrt(a2)*sqrt(b2)))
    
    return(alpha)
    
  }else{
    stop('Error, only two fp required ')
  } 
}

calcuate.fp.dist = function(xx0, use.cyst.radius, refine.rc = FALSE, use.global.rc = FALSE)
{
  # xx0 = res.cp[kk.fp, ]
  
  # first check if cyst id is unique
  if(length(unique(xx0$ID_cyst) == 1)){
    x0 = xx0[, grep('_cyst$', colnames(xx0))]
    yy = xx0[, grep('_fp$', colnames(xx0))]
    colnames(x0) = gsub('_cyst$', '', colnames(x0))
    colnames(yy) = gsub('_fp$', '', colnames(yy))
    extremeValue = 0
    
    if(use.global.rc){
      ## use the cyst radius to calculate the distance between fp
      dist.ratio = xx0$Distance_Centroid_organoid_fp/(xx0$AreaShape_EquivalentDiameter_cyst[1]/2.0)
      rc = xx0$Distance_Centroid_organoid_fp
      
      if(refine.rc){
        # filtering outliers of fp.to.cyst.distance with ratio thresholding if there are > 3 fp
        if(length(rc>=3)){
          rc = rc[which(dist.ratio>=0.5 & dist.ratio <=1.2)] 
        }
        
        if(length(rc)<1){
          extremeValue = 1
          #cat(xx0$ID_cyst[1],  ' fp.to.cyst.distance have extreme ratios\n')
          rc = xx0$Distance_Centroid_organoid_fp    
        }
      }
      
      rc = median(rc)
      #plot(res$AreaShape_Volume_fp.log10, res$Distance_Centroid_organoid_fp, cex = 0.2)
      #rrs = xx0$Distance_Centroid_organoid_fp
      #if(max(rrs)/min(rrs) > 2) cat('big difference in rc : ', unique(xx0$ID_cyst), '\n')
    }
   
    if(nrow(xx0) == 1){
      #rc = mean(calculate.distance(x0[1, ], yy))
      if(!use.global.rc) rc = xx0$Distance_Centroid_organoid_fp
      d.fp = 2*pi*rc
      return(c(d.fp, extremeValue))
      
    }else{
      
      if(nrow(xx0) == 2){
        if(!use.global.rc) rc = median(xx0$Distance_Centroid_organoid_fp)
        alpha = calculate.angle.between.two.fps(x0[1,], yy)
        d.fp = rc*alpha
        return(c(d.fp, extremeValue))
        
      }else{
        dds.fp = c()
        for(i in 1:nrow(yy))
        {
          dd.i = c()
          for(j in 1:nrow(yy))
          {
            if(j != i){
              if(!use.global.rc){
                dd.i = c(dd.i, mean(c(xx0$Distance_Centroid_organoid_fp[i], xx0$Distance_Centroid_organoid_fp[j]))*
                           calculate.angle.between.two.fps(x0[1,], yy[c(i, j), ]))
              }else{
                dd.i = c(dd.i, rc*calculate.angle.between.two.fps(x0[1,], yy[c(i, j), ])) 
              }
            }
          }
          
          dd.i = dd.i[order(dd.i)]
          dds.fp = c(dds.fp, mean(dd.i[c(1:2)]))
          
        }
        
        return(c(mean(dds.fp), extremeValue))
        
      } 
    }
    
    # rfp = (as.numeric(yy$Volume.Unit._Volume)/(4/3*pi))^(1/3)
    # rrc = median(dd + rfp)
    # return(rrc)
    
  }else{
    stop('multiple cyst ID found !!!')
  }
  
}

extract.turing.parameters.cellProfiler = function(res.cp, cyst.overlapRatio.threshold = 0.01, pixel.scale = 3,
                                                  refine.rc = FALSE)
{
  # res.cp = res
  cyst.id = unique(res.cp$ID_cyst)
  cc = res.cp$condition[match(cyst.id, res.cp$ID_cyst)]
  params =  data.frame(matrix(NA, nrow = length(cyst.id), ncol = 15))
  colnames(params) = c('nb.fp', 'dist.fp', 'dist.fp.extremeValue', 'foxa2.fp',  'radius.cyst', 'radius.fp', 'radiusSE.fp',
                       'volume', 'surface.area', 'overlap.ratio', 'olig2', 'foxa2', 'volume.fp', 'volumeSE.fp', 'dist.cyst.fp')
  rownames(params) = cyst.id
  
  # plot(res$AreaShape_EquivalentDiameter_fp/2, (res$Distance_Centroid_organoid_fp )/(res$AreaShape_EquivalentDiameter_cyst/2), 
  #      cex = 0.3, log = 'x');
  # #abline(h = c( 0.8), col ='red', lwd = 2.0)
  # rr = 10^seq(0, 5, length.out = 100)
  # points(rr, 1-rr/27, type = 'l', col = 'blue', lwd = 4.0)
  
  for(n in 1:length(cyst.id))
  {
    # n = 677; n = which(cyst.id == '107_11')
    if(n%%100 == 0) cat('cyst nb : ', n, '\n')
    kk = which(res.cp$ID_cyst == cyst.id[n])
    
    # extract information of cyst
    params$volume[n] = as.numeric(res.cp$AreaShape_Volume_cyst[kk[1]])
    params$surface.area[n] = as.numeric(res.cp$AreaShape_SurfaceArea_cyst[kk[1]])
    
    if(length(which(colnames(res) == 'Intensity_MeanIntensity_Olig2_cyst')) == 1 &
       length(which(colnames(res) == 'Scaling_Olig2_image')) == 1){
      params$olig2[n] = as.numeric(res.cp$Intensity_MeanIntensity_Olig2_cyst[kk[1]]) * res.cp$Scaling_Olig2_image[kk[1]]
    }
    
    params$foxa2[n] = as.numeric(res.cp$Intensity_MeanIntensity_FOXA2_cyst[kk[1]]) * res.cp$Scaling_FOXA2_image[kk[1]]
    params$radius.cyst[n] = as.numeric(res.cp$AreaShape_EquivalentDiameter_cyst[kk[1]])/2 # rough estimation with cyst volume
    
    kk.fp = kk[!is.na(res.cp$ID_fp[kk])]
    if(length(kk.fp) > 0){
      fp.cyst.ratio = sum(as.numeric(res.cp$AreaShape_Volume_fp[kk.fp]))/params$volume[n]
    }else{
      fp.cyst.ratio = 0
    }
    
    if(length(kk.fp) > 0 & fp.cyst.ratio > cyst.overlapRatio.threshold){
      params$nb.fp[n] = length(kk.fp)
      params$overlap.ratio[n] = sum(as.numeric(res.cp$AreaShape_Volume_fp[kk.fp]))/params$volume[n]
      
      #radius[n] = calcuate.cyst.radius(res.cp[kk.fp, ]) # calcuate this using the cyst center, fp center and fp radius
      params$foxa2.fp[n] = median(as.numeric(res.cp$Intensity_MeanIntensity_FOXA2_fp[kk.fp])) * res.cp$Scaling_FOXA2_image[kk[1]]
      
      params$volume.fp[n] = median(as.numeric(res.cp$AreaShape_Volume_fp[kk.fp]))
      params$radius.fp[n] = median(as.numeric(res.cp$AreaShape_EquivalentDiameter_fp[kk.fp]/2))
      params$dist.cyst.fp[n] = median(as.numeric(res.cp$Distance_Centroid_organoid_fp[kk.fp]))
      
      if(length(kk.fp)>1){
        params$volumeSE.fp[n] = sd(as.numeric(res.cp$AreaShape_Volume_fp[kk.fp]))
        params$radiusSE.fp[n] = sd(as.numeric(res.cp$AreaShape_EquivalentDiameter_fp[kk.fp]/2))
      }
      
      cyst.fp.dist = calcuate.fp.dist(res.cp[kk.fp, ], refine.rc = refine.rc)
      params$dist.fp[n] = cyst.fp.dist[1]
      params$dist.fp.extremeValue[n] = cyst.fp.dist[2]
      
    }else{
      params$nb.fp[n] = 0
      params$overlap.ratio[n] = sum(as.numeric(res.cp$AreaShape_Volume_fp[kk.fp]))/params$volume[n]
    }
  }
  
  # convert pixels to um
  params$dist.fp = params$dist.fp * pixel.scale
  params$radius.cyst = params$radius.cyst * pixel.scale
  params$radius.fp = params$radius.fp * pixel.scale
  params$dist.cyst.fp = params$dist.cyst.fp * pixel.scale
  params$radiusSE.fp = params$radiusSE.fp*pixel.scale
  params$volume = params$volume*pixel.scale^3
  params$volume.fp = params$volume.fp * pixel.scale^3
  params$volumeSE.fp = params$volumeSE.fp*pixel.scale^3
  
  
  params$surface.area = params$surface.area*pixel.scale^2/(4/3) # unknow factor from CP
  
  params = data.frame(condition = cc, params, stringsAsFactors = TRUE)
  
  for(n in 1:ncol(params))
  {
    if(colnames(params)[n] == 'condition'|colnames(params)[n] == 'nb.fp'){
      params[ ,n] = as.factor(params[,n])
    }else{
      params[,n] = as.numeric(params[,n])
    }
  }
  
 
  return(params)
  
}


########################################################
########################################################
# Section : merge tables from segmentation results by Imaris (NOT USED anymore !!!)
# 
########################################################
########################################################
make_mergedTables_fromSegementation_Imaris = function()
{
  cyst.channel = '1'
  floorplat.channel = '2'
  
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
  
  
}

extract.turing.parameters.imaris = function(res.cp)
{
  # res.cp = res
  cc = unique(res.cp$condition)
  params = c()
  
  for(m in 1:length(cc))
  {
    # m = 18
    c = cc[m]
    cat('condition : ', c, '\n')
    
    # parameters to extract
    cyst.id = unique(res.cp$ID_cyst[which(res.cp$condition == c)])
    
    volume = rep(NA, length(cyst.id))
    area = rep(NA, length(cyst.id))
    voxel = rep(NA, length(cyst.id))
    overlap.ratio = rep(NA, length(cyst.id))
    olig2 = rep(NA, length(cyst.id))
    foxa2 = rep(NA, length(cyst.id))
    
    radius = rep(NA, length(cyst.id))
    nb.fp = rep(NA, length(cyst.id))
    
    dist.fp = rep(NA, length(cyst.id))
    volume.fp = rep(NA, length(cyst.id))
    radius.fp = rep(NA, length(cyst.id))
    foxa2.fp = rep(NA, length(cyst.id))
    
    for(n in 1:length(cyst.id))
    {
      # n = 4
      kk = which(res.cp$ID_cyst == cyst.id[n] & res.cp$condition == c)
      
      # extract information of cyst
      volume[n] = as.numeric(res.cp$Volume.Unit._Volume_cyst[kk[1]])
      area[n] = as.numeric(res.cp$Area.Unit.µm2_Area_cyst[kk[1]])
      voxel[n] = as.numeric(res.cp$Number.of.Voxels.Unit._Number.of.Voxels.Img1_cyst[kk[1]])
      overlap.ratio[n] = as.numeric(res.cp$Overlapped.Volume.Ratio_cyst[kk[1]])
      olig2[n] = as.numeric(res.cp$Intensity.Mean.Unit._Intensity.Mean.Ch2.Img1_cyst[kk[1]])
      foxa2[n] = as.numeric(res.cp$Intensity.Mean.Unit._Intensity.Mean.Ch3.Img1_cyst[kk[1]])
      
      kk.fp = kk[!is.na(res.cp$ID_fp[kk])]
      #cat(length(kk), 'fp \n')
      
      if(length(kk.fp) == 0){ # no fp in the cyst
        nb.fp[n] = 0
        radius[n] = ((volume[n])/(4/3*pi))^(1/3) # rough estimation with cyst volume
      }else{
        nb.fp[n] = length(kk.fp)
        radius[n] = calcuate.cyst.radius(res.cp[kk.fp, ]) # calcuate this using the cyst center, fp center and fp radius
        volume.fp[n] = median(as.numeric(res.cp$Volume.Unit._Volume_fp[kk.fp]))
        foxa2.fp[n] = median(as.numeric(res.cp$Intensity.Mean.Unit._Intensity.Mean.Ch3.Img1_fp[kk.fp]))
        radius.fp[n] = median((volume.fp[n])/(4/3*pi)^(1/3))
        
        dist.fp[n] = calcuate.fp.dist(res.cp[kk.fp, ])
      }
      
    }
    
    params = rbind(params, cbind(rep(c, length(cyst.id)), cyst.id, nb.fp, dist.fp, foxa2.fp, radius.fp, radius, 
                                 volume, area, voxel, overlap.ratio,
                                 olig2, foxa2, volume.fp))
    
  }
  
  params = data.frame(params, stringsAsFactors = FALSE)
  colnames(params) = c('condition', 'cyst.id', 'nb.fp', 'dist.fp', 'foxa2.fp', 'radius.fp',  'radius.cyst', 
                       'volume', 'area', 'voxel', 'overlap.ratio', 'olig2', 'foxa2', 'volume.fp')
  
  return(params)
  
}




find.metadata.from.imaris = function()
{
  xx = readRDS(file = paste0('../results/Rdata/RdatamergedTable_cyst.fp_allConditions.rds'))
  design = data.frame(xx$condition, xx$Original.Image.Name_cyst, xx$Original.Image.ID_cyst, stringsAsFactors = FALSE)
  colnames(design) = c('condition', 'Original.Image.Name', 'Original.Image.ID')
  image.uniq = unique(design$Original.Image.Name)
  
  design = design[match(image.uniq, design$Original.Image.Name), ]
  
  saveRDS(design, file = paste0(Rdata, '/perturbation_design_hNTdrugs3_0310.rds'))
  
  
}

test.CPversion.output.shape.metric = function()
{
  # test the cyst from human organoid data
  cyst = read.csv(file = paste0(dataDir, 'organoid.csv'))
  xx = cyst
  
  xx$sphericity = pi^(1/3)*(6*xx$AreaShape_Volume)^(2/3)/xx$AreaShape_SurfaceArea
  xx$R0 = sqrt(xx$AreaShape_SurfaceArea/(4*pi))
  xx$reduced.volume = xx$AreaShape_Volume/(4/3*pi*xx$R0^3)
  
  par(mfrow=c(1,2)) 
  plot(xx$AreaShape_Volume, xx$sphericity, cex = 0.3, ylim = c(0, 1), main = 'volume and sphericity',
       xlab = 'AreaShape_Volume', ylab = 'Sphericity')
  abline(h = 3/4, col = 'red', lwd = 2.0)
  
  plot(xx$AreaShape_Volume, xx$reduced.volume, cex = 0.3, ylim = c(0, 1), main = 'volume and reduced volume', 
       xlab = 'AreaShape_Volume', ylab = 'reduced.volume')
  abline(h = (3/4)^(3/2), col = 'red', lwd = 2.0)
  
  write.table(xx, file = paste0(tabDir, '/cyst_shapeTest_sphericity_reducedVolume.txt'), sep = '\t', col.names = TRUE, row.names = FALSE)
  
  
  ##########################################
  # test surface and volume in CP with perfect spheres 
  ##########################################
  image = read.csv(file = '/Users/jiwang/workspace/imp/organoid_patterning/images/cellProfiler_test/sphere_test2/MyExpt_Image.csv')
  xx = read.csv(file = '/Users/jiwang/workspace/imp/organoid_patterning/images/cellProfiler_test/sphere_test2/MyExpt_organoid.csv')
  xx = data.frame(image[match(xx$ImageNumber, image$ImageNumber)], xx)
  xx$volume =  xx$AreaShape_Volume
  xx$surface = xx$AreaShape_SurfaceArea
  xx$radius = xx$AreaShape_EquivalentDiameter/2
  #xx$sphericity = pi^(1/3)*(6*xx$volume)^(2/3)/xx$surface
  rr = c(100, 10, 200, 25, 50)
  
  par(mfrow=c(1,3)) 
  plot(rr, xx$radius, xlab = 'true values', ylab = 'estimated from CP', main = 'Rdadius', cex = 1.2); 
  abline(0, 1, lwd = 2.0, col = 'red')
  plot(4/3*pi*rr^3, xx$volume, xlab = 'true values', ylab = 'estimated from CP', main = 'Volume', log = 'xy', cex = 1.2);
  abline(0, 1, lwd = 2.0, col = 'red')
  plot(4*pi*rr^2, xx$surface, xlab = 'true values', ylab = 'estimated from CP', main = 'Surface area', log = 'xy', cex = 1.2);
  points(seq(1, 10^6, 200), seq(1,10^6, 200)*4/3, type = 'l', col = 'blue')
  abline(0, 1, lwd = 2.0, col = 'red')
  #abline(0, 4/3, lwd = 2.0, col = 'red', col = )
  
  
  
  plot(xx$volume, xx$sphericity)
  
  plot(xx$Volume..pixel.3.)
  
  plot(4/3*4*pi*(xx$AreaShape_EquivalentDiameter/2)^2, xx$AreaShape_SurfaceArea, cex = 0.6)
  abline(0, 1, col = 'red')
  
  
}

##########################################
# Filtering cyst and fp after concatenating the table 
##########################################
doubleCheck.CP.surfaceArea.volume = function(xx)
{
  plot(xx$AreaShape_Volume_cyst, xx$AreaShape_SurfaceArea_cyst)
  rr = c(0:1000)
  points(4/3*pi*rr^3, 4*pi*rr^2, type = 'l')
  points(rr^3, 6*rr^2, type = 'l')
  
  plot(4*pi*(xx$AreaShape_EquivalentDiameter_cyst/2)^2, xx$AreaShape_SurfaceArea_cyst, cex = 0.6)
  abline(0, 1, col = 'red')
  
  plot(4/3*pi*(xx$AreaShape_EquivalentDiameter_cyst/2)^3, xx$AreaShape_Volume_cyst);
  abline(0, 1, lwd =2.0, col = 'red')
}





