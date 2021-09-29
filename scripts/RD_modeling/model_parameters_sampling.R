##########################################################################
##########################################################################
# Project: NT organoid patterning 
# Script purpose: genereate the network topologies for RD 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Sep  1 16:37:34 2021
##########################################################################
##########################################################################
rm(list = ls())

library(igraph)

network = '4N3M'


if(network == '3N2M'){
  resDir = '../results/RD_topology_screening/3N2M_topology_enumerate_v1/'
  if(!dir.exists(resDir)) dir.create(resDir)
  
  s0 = matrix(NA, nrow = 3, ncol = 3)
  colnames(s0) = c('Nog', 'BMP', 'Foxa2')
  rownames(s0) = colnames(s0)
  
  s0[2, 1] = 1 # BMP activate Nog from TM3 data
  s0[3, 1] = 1 # Foxa2 cells expressing Nog
  s0[1, 2] = -1; # Nog inhibite BMP
  s0[3, 3] = 1 # Foxa2 self activation
  
  
  ii2assign = which(is.na(s0))
  
  xx = expand.grid(0:1, # noggin auto-activation 
                   -1:1, 
                   -1:1, 
                   0:1, # Nog activate Foxa2
                   -1:1)
  
}

if(network == '4N3M')
{
  resDir = 'RD_modeling/4N3M_topology_enumerate/'
  if(!dir.exists(resDir)) dir.create(resDir)
  
  s0 = matrix(NA, nrow = 4, ncol = 4)
  colnames(s0) = c('Nog', 'BMP', 'Shh', 'Foxa2')
  rownames(s0) = colnames(s0)
  s0[1, 1] = 0 # no Noggin autoactivation
  s0[1, 2] = -1; # Nog inhibite BMP
  s0[1, 3] = 0; # no interaction between Shh and Noggin
  s0[1, 4] = 0; # no positive feedback from noggin to foxa2
  
  s0[2, 1] = 1 # BMP activate Nog from TM3 data
  s0[2, 4] = -1 # bmp inhibites foxa2
  
  s0[3, 1] = 0 # no interaction between Shh and Noggin
  s0[3, 3] = 0 # no autoregulation of Shh
  s0[3, 4] = 1 # shh activate foxa2
  
  s0[4, 1] = 1 # Foxa2 activate noggin
  s0[4, 3] = 1# foxa2 activate shh
  
  
  ii2assign = which(is.na(s0))
  
  xx = expand.grid(-1:1, # bmp autoregulation
                   -1:1, # shh to bmp
                   0:1, # Foxa2 expresses or not bmp
                   -1:1, # bmp regultes shh
                   0:1)  # foxa2 auto-activation
  
}


#nb_model = 1
for(n in 1:nrow(xx))
{
  # n = 1
  cat(n, '\n')
  s = s0
  s[ii2assign] = unlist(xx[n, ])
  
  write.csv(s, file = paste0(resDir, 'Model_', n, '.csv'), row.names = TRUE, quote = FALSE)
  
  pdfname = paste0(resDir, 'Model_', n,  ".pdf")
  pdf(pdfname, width = 8, height = 6)
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  
  make.plot.from.adjacency.matrix(as.matrix(s), main = paste0('Model_', n), network = network)
  
  dev.off()
  
  #g1 = graph_from_adjacency_matrix(s, mode = 'directed')
  #plot(g1)
  
}

########################################################
########################################################
# Section : summary of network screening
# 
########################################################
########################################################
rm(list = ls())
library("igraph")

RDoutDir = '../results/RD_topology_screening/topology_screening_4N3M_v1/'
screening.outDir = paste0(RDoutDir, 'RD_out_4N3M/')
modelDir = paste0(RDoutDir, '4N3M_topology_enumerate')

# make graph from adjacency matrix
make.plot.from.adjacency.matrix = function(s, main = '', network = '3N2M')
{
  g1 = graph_from_adjacency_matrix(s, mode = 'directed', weighted =TRUE)
  ws = E(g1)$weight
  cols = rep('darkblue', length(ws))
  cols[which(ws<0)] = 'red'
  if(network == '3N2M') {
    layout = matrix(data = c(0, 4, 2, 0, 0, 2), ncol = 2)
  }else{
    layout = matrix(data = c(0, 4, 4, 0, 0, 0, 2, 2), ncol = 2)
  }
  
  plot(g1, layout = layout,  layout = layout.circle, vertex.size = 60.0, edge.curved=.3, edge.color=cols, edge.arrow.size=2.0,
       edge.width = 4.0, vertex.label.cex = 1.8, vertex.label.font = 2, main = main)
  
}

filter.reaction.parameters.for.RD.patterning = function(res, filter.phase = FALSE, filter.oscillation = FALSE, network = '3N2M')
{
  lambda_real = res[ , grep('lambda_re', colnames(res))]
  lambda_im = res[ , grep('lambda_im', colnames(res))]
  k = res[, match(paste0('q', c(0:(ncol(lambda_real)-1))) , colnames(res))]
  phases =  res[ , grep('eigenvec_posSign', colnames(res))]
  
  index_keep = c()
  for(jj in 1:nrow(res))
  {
    # remove the turing filtering by checking if the Re(lambda) is negative for wavenumber max
    # BMP diffusion is >=0.1, comparable or larger than Noggin
    if(!filter.phase & !filter.oscillation){
      if(lambda_real[jj, which(k[jj, ] == max(k[jj, ]))] < 0  &  res$d1[jj] >= 0.5) index_keep = c(index_keep, jj)
    }
    if(filter.phase){
      pp = phases[jj, which(lambda_real[jj, ]>0)]
      
      if(network == '3N2M'){
        if(lambda_real[jj, which(k[jj, ] == max(k[jj, ]))] < 0  &  res$d1[jj] >= 0.5 & (all(pp == '1;0;1') | all(pp == '0;1;0'))) {
          index_keep = c(index_keep, jj)
        }
      }
      if(network == '4N3M') {
        if(lambda_real[jj, which(k[jj, ] == max(k[jj, ]))] < 0  &  res$d1[jj] >= 0.5 & (all(pp == '1;0;1;1'| pp == '0;1;0;0'))) {
          index_keep = c(index_keep, jj)
        }
      }
      
    }
    
  }
  
  return(res[index_keep, ])
  
}

s0 = matrix(NA, nrow = 3, ncol = 3)
s0[1, 1] = 0; s0[1, 2] = -1; s0[1, 3] = 0;
s0[2, 1] = 1; s0[2, 2] = 0; s0[2, 3] = -1
s0[3, 1] = 1; s0[3, 2] = 1; s0[3, 3] = 1

colnames(s0) = c('Nog', 'BMP', 'Foxa2')
rownames(s0) = colnames(s0) 
# make.plot.from.adjacency.matrix(s0)

model.list = list.files(path = modelDir, pattern = '*.csv', full.names = TRUE)
network = '4N3M'

nb = 0

#selection.criterion = 'D.larger.1_lambda.neg.max.q_phase'
resDir = paste0(RDoutDir, 'topology_summary_selection_D.larger.1_lambda.neg.max.q_phase')
modelSaveDir = paste0(resDir, '/Models_selected')
paramSaveDir = paste0(resDir, '/table_params')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(modelSaveDir)) dir.create(modelSaveDir)
if(!dir.exists(paramSaveDir)) dir.create(paramSaveDir)

for(n in 1:length(model.list)){
  # n = 94
  Model = basename(model.list[n])
  Model = gsub('.csv', '', Model)
  
  param.list = list.files(path = paste0(screening.outDir, Model), 
                          pattern = '*.csv', full.names = TRUE)
  
  if(length(param.list) > 0){
    
    ss = read.csv(model.list[n], header = TRUE, row.names = 1)
    
    nb.param = 0
    index.param = c()
    keep = c()
    for(i in 1:length(param.list))
    {
      # i = 3
      res = read.csv(file = param.list[i], header = TRUE)
      
      ##########################################
      # 1st filter the parameters for which Re(lambda) > 0 when q = 0
      ##########################################
      res = res[which(res$noDiffusion0 <0 ), ]
      res = filter.reaction.parameters.for.RD.patterning(res, filter.phase = TRUE, network = network)
      
      if(nrow(res) > 0) {
        nb.param = nb.param + 1
        index.param = c(index.param, rep(nb.param, nrow(res)))
        keep = rbind(keep, res)
      }
    }
    
    # cat(Model, '---', nb.param, 'sets of reaction parameters\n')
    
    keep = data.frame(index.param = index.param, keep, stringsAsFactors = FALSE)
    
    if(nb.param >0) {
      nb = nb + 1
      cat('Model index: ',  n, ': ', Model,  ' -- nb of reaction params:',  nb.param, '\n')
      
      pdfname = paste0(modelSaveDir, "/", Model, ".pdf")
      pdf(pdfname, width = 8, height = 6)
      par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
      
      if(network == '3N2M'){
        make.plot.from.adjacency.matrix(as.matrix(ss), 
                                        main = paste0(Model, ' , Q = ', nb.param, ', D_min = ', signif(min(keep$d1), d=2)),
                                        network = network)
      }
      if(network == '4N3M'){
        make.plot.from.adjacency.matrix(as.matrix(ss), 
                                        main = paste0(Model, ' , Q = ', nb.param, ', D_bmp.shh.min = ', signif(min(keep$d1), d=2), 
                                                      signif(min(keep$d2), d=2)),
                                        network = network)
      }
      
      dev.off()
      
      # make summary for the model
      pdfname = paste0(resDir, "/", Model, ".pdf")
      pdf(pdfname, width = 12, height = 8)
      par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
      
      #s = as.matrix(ss)
      
      if(network == '3N2M'){
        make.plot.from.adjacency.matrix(as.matrix(ss), 
                                        main = paste0(Model, ' , Q = ', nb.param, ', D_min = ', signif(min(keep$d1), d=2)),
                                        network = network)
      }
      if(network == '4N3M'){
        make.plot.from.adjacency.matrix(as.matrix(ss), 
                                        main = paste0(Model, ' , Q = ', nb.param, ', D_bmp.shh.min = ', signif(min(keep$d1), d=2), 
                                                      signif(min(keep$d2), d=2)),
                                        network = network)
      }
      
      for(j in unique(keep$index.param))
      {
        # j = 1
        sels = which(keep$index.param == j)
        
        lambda_real = keep[sels, grep('lambda_re', colnames(keep))]
        lambda_im = keep[sels, grep('lambda_im', colnames(keep))]
        k = keep[sels, match(paste0('q', c(0:(ncol(lambda_real)-1))) , colnames(keep))]
        
        plot(1, 1, type = 'n', xlim = range(k), ylim = range(lambda_real), log='x', xlab = 'wavenumber (k)', ylab = 'Re(lamba.max)',
             main = Model)
        abline(h = 0, lwd = 2.0, col = 'gray80')
        for(ii in 1:nrow(k)){
          points(as.numeric(k[ii, ]), as.numeric(lambda_real[ii, ]), 
               type = 'l', col = ii, lwd = 2.0)
          
        }
        if(network == '3N2M'){
          legend('topleft', legend = paste0('d = ', signif(keep$d1[sels], d = 2)), 
                 col = 1:nrow(k), cex = 0.7)
        }
        if(network == '4N3M'){
          legend('topleft', legend = paste0('d = ', signif(keep$d1[sels], d = 2), ';', signif(keep$d2[sels], d = 2)), 
                 col = 1:nrow(k), cex = 0.7)
        }
        
      }
      
      dev.off()
      
      write.csv(keep, file = paste0(paramSaveDir, '/params_saved_', Model, '.csv'), row.names = FALSE, quote = FALSE)
      write.csv(ss, file = paste0(paramSaveDir, '/', Model, '.csv'), row.names = TRUE, quote = FALSE)
      
    }  
  }
}
