##########################################################################
##########################################################################
# Project: NT organoid patterning 
# Script purpose: genereate the network parameters for RD 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Sep  1 16:37:34 2021
##########################################################################
##########################################################################
resDir = paste0("RD_modeling/3N2M_topology_enumerate/")
if(!dir.exists(resDir)) dir.create(resDir)

s = matrix(NA, nrow = 3, ncol = 3)
colnames(s) = c('Nog', 'BMP', 'Foxa2')
rownames(s) = colnames(s)
s[1, 1] = 0; s[1, 2] = -1; s[1, 3] = 0;
s[2, 1] = 1; s[2, 2] = 0; s[2, 3] = -1
s[3, 1] = 1; s[3, 2] = 1; s[3, 3] = 1

write.csv(s, file = paste0(resDir, 'Model_1.csv'), row.names = TRUE, quote = FALSE)

# library(emdbook)
# library(tidyr)
# nb_sampling_perParam = 3
# ks = lseq(0.1, 100, length.out = 3)
# gamma = lseq(0.01, 1, length.out = 3)
# 
# tidyr::expand(data.frame(ks), data.frame(gamma))
