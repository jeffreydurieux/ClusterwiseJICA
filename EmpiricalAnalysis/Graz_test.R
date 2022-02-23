# Wed Feb 23 10:07:51 2022
# Author: Jeffrey Durieux, MSc

# What: script empirical test on graz data

library(CICA)
library(mclust)
library(plotly)

source('sortX.R')
source('ICAonList.R')
source('computeAhats.R')
source('computeXhats.R')
source('Avoid_nc_N.R')
source('Simulate_CJICA.R')
source('ClusterwiseJICA.R')


load('~/Documents/Grazdata_CJICA/AlzheimerfMRI_extra_aangevuld.RData')

#### Take VBM4mm and ALFF4mm --> Christiane's best option

dim(ALFF4mm)  # 25750
dim(VBM4mm)   # 59049

X <- cbind(VBM4mm, ALFF4mm) # 250 x 59049 | 25750



