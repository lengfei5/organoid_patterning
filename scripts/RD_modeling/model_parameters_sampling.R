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
resDir = paste0("RD_modeling/3N2M_topology_enumerate/")
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
nb_model = 1

for(n in 1:nrow(xx))
{
  # n = 21
  cat(n, '\n')
  s = s0
  s[ii2assign] = unlist(xx[n, ])
  
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
library("igraph")

s0 = matrix(NA, nrow = 3, ncol = 3)
s0[1, 1] = 0; s0[1, 2] = -1; s0[1, 3] = 0;
s0[2, 1] = 1; s0[2, 2] = 0; s0[2, 3] = -1
s0[3, 1] = 1; s0[3, 2] = 1; s0[3, 3] = 1

colnames(s0) = c('Nog', 'BMP', 'Foxa2')
rownames(s0) = colnames(s0)

# make graph from adjacency matrix
g1 = graph_from_adjacency_matrix(s0, mode = 'directed', weighted =TRUE)
edge.start <- ends(g1, es=E(g1), names=F)[,1]
edge.col <- V(g1)$color[edge.start]
plot(g1, edge.color=edge.col, edge.curved=.1)  
plot(g1, edge.color = 'blue')



L3 <- LETTERS[1:8]
d <- data.frame(start = sample(L3, 16, replace = T), end = sample(L3, 16, replace = T),
                weight = c(20,40,20,30,50,60,20,30,20,40,20,30,50,60,20,30))


g <- graph.data.frame(d, directed = T)
V(g)$name 
E(g)$weight

ideg <- degree(g, mode = "in", loops = F)

col=rainbow(12) # For edge colors

plot.igraph(g, 
            vertex.label = V(g)$name, vertex.label.color = "gray20",
            vertex.size = ideg*25 + 40, vertex.size2 = 30,
            vertex.color = "gray90", vertex.frame.color = "gray20",
            vertex.shape = "rectangle",
            edge.arrow.size=0.5, edge.color=col, edge.width = E(g)$weight / 10,
            edge.curved = T, 
            layout = layout.reingold.tilford)

gridLayout <- function(x)
{
  LmatX <- seq(-1,1,length=ncol(x))
  LmatY <- seq(1,-1,length=nrow(x))
  
  loc <- t(sapply(1:max(x),function(y)which(x==y,arr.ind=T)))
  layout <- cbind(LmatX[loc[,2]],LmatY[loc[,1]])
  return(layout)
}

grid <- matrix(c(
  0,0,1,0,0,
  2,0,3,0,4),nrow=2,byrow=TRUE)

library("igraph")

g <- graph.adjacency(matrix(1,4,4))

plot(g,layout=gridLayout(L))


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


