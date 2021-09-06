##########################################################################
##########################################################################
# Project: NT organoid patterning 
# Script purpose: genereate the network parameters for RD 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Sep  1 16:37:34 2021
##########################################################################
##########################################################################
library(emdbook)
library(tidyr)
nb_sampling_perParam = 3
ks = lseq(0.1, 100, length.out = 3)
gamma = lseq(0.01, 1, length.out = 3)

tidyr::expand(data.frame(ks), data.frame(gamma))