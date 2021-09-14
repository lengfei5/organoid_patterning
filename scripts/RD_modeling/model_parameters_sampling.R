##########################################################################
##########################################################################
# Project: NT organoid patterning 
# Script purpose: genereate the network topologies for RD 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Sep  1 16:37:34 2021
##########################################################################
##########################################################################
library(igraph)
resDir = paste0("RD_modeling/3N2M_topology_enumerate/")
if(!dir.exists(resDir)) dir.create(resDir)

s = matrix(NA, nrow = 3, ncol = 3)
colnames(s) = c('Nog', 'BMP', 'Foxa2')
rownames(s) = colnames(s)

s[2, 1] = 1 # BMP activate Nog from TM3 data
s[3, 1] = 1 # Foxa2 cells expressing Nog
s[1, 2] = -1; # Nog inhibite BMP
s[3, 3] = 1 # Foxa2 self activation


ii.toassing = which(is.na(s))

xx = expand.grid(0:1, 0:1, -1:1, -1:1, -1:1)
nb_model = 1

for(n in 1:nrow(xx))
{
  # n = 1
  cat(n, '\n')
  s[ii.toassing] = unlist(xx[n, ])
  
  write.csv(s, file = paste0(resDir, 'Model_', n, '.csv'), row.names = TRUE, quote = FALSE)
  
  #g1 = graph_from_adjacency_matrix(s, mode = 'directed')
  #plot(g1)
  
}

########################################################
########################################################
# Section : summary of models 
# 
########################################################
########################################################
s = matrix(NA, nrow = 3, ncol = 3)
s[1, 1] = 0; s[1, 2] = -1; s[1, 3] = 0;
s[2, 1] = 1; s[2, 2] = 0; s[2, 3] = -1
s[3, 1] = 1; s[3, 2] = 1; s[3, 3] = 1

s0 = s
# make graph from adjacency matrix
g1 = graph_from_adjacency_matrix(s, mode = 'directed', weighted =TRUE)
edge.start <- ends(g1, es=E(g1), names=F)[,1]
edge.col <- V(g1)$color[edge.start]
plot(g1, edge.color=edge.col, edge.curved=.1)  
plot(g1, edge.color = 'blue')


# library(emdbook)
# library(tidyr)
# nb_sampling_perParam = 3
# ks = lseq(0.1, 100, length.out = 3)
# gamma = lseq(0.01, 1, length.out = 3)
# 
# tidyr::expand(data.frame(ks), data.frame(gamma))

resDir = paste0("../results/RD_topology_test/topology_screening_3N2M")
if(!dir.exists(resDir)) dir.create(resDir)


model.list = list.files(path = "RD_modeling/3N2M_topology_enumerate", pattern = '*.csv', full.names = TRUE)
nb = 0

for(n in 1:length(model.list)){
  # n = 20
  Model = basename(model.list[n])
  Model = gsub('.csv', '', Model)
  
  param.list = list.files(path = paste0('../results/RD_topology_test/RD_out_topology.108_params.50K/', Model), 
                          pattern = '*.csv', full.names = TRUE)
  
  if(length(param.list) > 0){
    
    
    ss = read.csv(model.list[n], header = TRUE, row.names = 1)
    nb.param = 0
    index.param = c()
    keep = c()
    for(i in 1:length(param.list))
    {
      # i = 15
      res = read.csv(file = param.list[i], header = TRUE)
      res = res[which(res$noDiffusion0 <0 ), ]
      if(nrow(res) > 0) {
        nb.param = nb.param + 1
        index.param = c(index.param, rep(nb.param, nrow(res)))
        keep = rbind(keep, res)
        cat(nb.param, 'th parameters \n')
      }
    }
    
    keep = data.frame(index.param = index.param, keep, stringsAsFactors = FALSE)
    
    if(nb.param >0) {
      nb = nb + 1
      cat(n, ' -- ', Model,  ': nb of param --', nb.param,  '-- nb : ', nb, '\n')
      
      # make summary for the model
      pdfname = paste0(resDir, "/", Model, ".pdf")
      pdf(pdfname, width = 12, height = 8)
      par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
      
      for(j in unique(keep$index.param))
      {
        # j = 1
        sels = which(keep$index.param == j)
        k = keep[sels, match(paste0('q', c(0:19)) , colnames(keep))]
        lambda_real = keep[sels, grep('lambda_re', colnames(keep))]
        lambda_im = keep[sels, grep('lambda_im', colnames(keep))]
        
        plot(1, 1, type = 'n', xlim = range(k), ylim = range(lambda_real), log='x', xlab = 'wavenumber (k)', ylab = 'Re(lamba.max)',
             main = Model)
        abline(h = 0, lwd = 2.0, col = 'gray80')
        for(ii in 1:nrow(k)){
          points(as.numeric(k[ii, ]), as.numeric(lambda_real[ii, ]), 
               type = 'l', col = ii, lwd = 2.0)
          
        }
        
        legend('topleft', legend = paste0('d = ', signif(keep$d1[sels], d = 2), col = 1:nrow(k)))
        
        
      }
      
      dev.off() 
      
    }  
  }
  
}


