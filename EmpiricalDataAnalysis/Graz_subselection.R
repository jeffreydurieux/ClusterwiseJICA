# Tue Jun  8 09:38:01 2021
# Author: Jeffrey Durieux, MSc

# test on Graz data 20 vs 20


load('~/Documents/Grazdata_CJICA/names_of_data.Rdata')
load('~/Documents/Grazdata_CJICA/subselectionnames.Rdata')

files <- gsub(pattern = '_.*', replacement = '', x = files)

id_sub <- names %in% files

load('~/Documents/Grazdata_CJICA/VBM4mm.Rdata')
load('~/Documents/Grazdata_CJICA/ALFF4mm.Rdata')

dim(ALFF4mm)
dim(VBM4mm)

ALFF <- ALFF4mm[id_sub, ]
VBM <- VBM4mm[id_sub, ]

# check ss and scale
sum(ALFF^2)
sum(VBM^2)

f1 <- sqrt(5000/sum(VBM^2))
f2 <- sqrt(5000/sum(ALFF^2))
X1 <- f1*VBM
X2 <- f2*ALFF
sum(X1^2)
sum(X2^2)

X <- cbind(X1,X2)
dim(X)
X <- t(X)

##### analysis ######
library(mclust) #ARI
library(plotly) # '>' operator
source('sortX.R')
source('ICAonList.R')
source('computeAhats.R')
source('computeXhats.R')
source('Avoid_nc_N.R')
source('Simulate_CJICA.R')
source('ClusterwiseJICA.R')
source('SearchEmptyClusters.R')
source('ILS_CJICA.R')
source('rankperb.R')
source('ssranking.R')
source('scaleprob_uij2.R')

lab <- c(rep(1,20), rep(2,20))
Q = 5
Grazsub_q2_r2 <- ClusterwiseJICA(X = X, k = 2,
                nc = Q, starts = 50, scale = F)


losses <- sapply(seq_along(Grazsub_q2_r2), 
       function(i) tail(Grazsub_q2_r2[[i]]$lossiter, n = 1))
plot(losses)
optid <- which.min(losses)

opt <- Grazsub_q2_r2[[optid]]
opt$lossiter
opt$Lir$loss
tab <- table(lab,opt$p)
tab
caret::confusionMatrix(tab)

ils <- ILSclusterwise(X = X, p = opt$p, R = 2, Q = Q)
ils$Lir$loss
ils$p
ils$loss_of_p

tab <- table(rev(lab),ils$p)
tab
caret::confusionMatrix(tab)
mclust::adjustedRandIndex(lab,ils$p)

ilstrue <- ILSclusterwise(X = X, p = lab, R = 2, Q = Q)
ilstrue$Lir$loss
ilstrue$p
ilstrue$loss_of_p
min(ilstrue$lossiter)

lab2 <- ifelse(lab==1,2,1)
tab <- table(rev(lab),ilstrue$p)
tab
caret::confusionMatrix(tab)
mclust::adjustedRandIndex(lab,ilstrue$p)

