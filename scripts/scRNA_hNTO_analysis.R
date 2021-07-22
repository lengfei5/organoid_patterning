##########################################################################
##########################################################################
# Project: organoid patterning project
# Script purpose: analyze the scRNA-seq data from Adrian Ranga https://www.nature.com/articles/s41467-021-22952-0
# Here the Day 3 (before RA and SAG) and Day 5 (right after RA washing and without RA and SAG) data were selected;
# Day 11 data is not being considered, becasue the pattern at day 11 is not clear as salt-and-pepper, because it does not seem to be 
# stable pattern. 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Thu Jul 22 13:50:47 2021
##########################################################################
##########################################################################

rm(list = ls())
RNA.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_functions.R'
RNA.QC.functions = '/Volumes/groups/tanaka/People/current/jiwang/scripts/functions/RNAseq_QCs.R'
source(RNA.functions)
source(RNA.QC.functions)

# setup for data import and sequencing QCs
version.analysis = '_v20210722'

resDir = paste0("../results/scRNAseq_hNTO_Adrian", version.analysis)
RdataDir = paste0(resDir, '/Rdata/')

if(!dir.exists(resDir)) dir.create(resDir)
if(!dir.exists(RdataDir)) dir.create(RdataDir)

dataDir = '/Volumes/groups/tanaka/People/current/jiwang/projects/patterning_organoid/R11601/'


