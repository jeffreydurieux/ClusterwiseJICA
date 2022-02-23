# Mon Oct 26 14:28:14 2020
# Author: Jeffrey Durieux, MSc
#

# some tests

source('sortX.R')
source('ICAonList.R')
source('computeAhats.R')
source('computeXhats.R')
source('Avoid_nc_N.R')
source('Simulate_CJICA.R')
source('ClusterwiseJICA.R')

library(CICA)
library(mclust)
library(plotly)
rm(test1)
#### Graz####
load('~/Documents/Grazdata_CJICA/AlzheimerfMRI_extra_aangevuld.RData')
X <- cbind(VBM4mm, FCdefaultModeNetwork4mm)
X <- t(X)
dim(X)
rm(correlations20components, correlations70components, partialCorrelations20components, partialCorrelations70components, FCdefaultModeNetwork2mm, FCdefaultModeNetwork4mm, FCexecutiveControl2mm, FCexecutiveControl4mm, VBM2mm, VBM4mm, ALFF2mm, ALFF4mm)

Alzheimer <- ifelse(Alzheimer == 0, 1, 2)
set.seed(2407)
graztest <- ClusterwiseJICA(X = X, k = 2, nc = 5, scale = F, starts = 1, rational = Alzheimer)
graztest[[1]]$p


table(Alzheimer)

table(graztest[[1]]$p, Alzheimer)
tab <- table(graztest[[1]]$p, Alzheimer)
caret::confusionMatrix(tab)


tt <- table(graztest[[1]]$p, Alzheimer)
caret::confusionMatrix(tt)


#### simulated data ######

K <- 2; Q <- 2; Nk = 50
data <- Simulate_CJICA(Nk = Nk, Vm = 2500, K = K, Qm = Q, E = .05
                       , M = 2, type = 4, cor = .71)
str(data)
cor(data$S[[1]],data$S[[2]])
#cor(data$S[[1]],data$S[[3]])
#modRV(data$S[[1]],data$S[[2]])

time <- proc.time()
cjica <- ClusterwiseJICA(X = data$Xe, k = K, nc = Q, starts = 10)
tm <- proc.time() - time 
tm

minloss <- sapply(seq_along(cjica), function(anom) tail(cjica[[anom]]$lossiter,1))
plot(minloss)

id <- minloss == min(minloss)

resmin <- cjica[id]

sols <- length(resmin)
id <- sample(sols,1)
rr <- resmin[[id]]
rr$p
tail(rr$lossiter,1)


ll <- tail(rr$lossiter,1)


adjustedRandIndex(data$P,rr$p)
cor(data$S[[1]], rr$ica$Sr[[1]]) %>% round(digits = 2)
cor(data$S[[2]], rr$ica$Sr[[2]]) %>% round(digits = 2)
cor(data$S[[1]], rr$ica$Sr[[2]]) %>% round(digits = 2)
cor(data$S[[2]], rr$ica$Sr[[1]]) %>% round(digits = 2)

modRV(data$S[[1]], rr$ica$Sr[[1]])
modRV(data$S[[2]], rr$ica$Sr[[2]])
modRV(data$S[[1]], rr$ica$Sr[[2]])
modRV(data$S[[2]], rr$ica$Sr[[1]])

cor(data$A[[1]], rr$ica$Mr[[1]]) %>% round(digits = 3)
cor(data$A[[1]], rr$ica$Mr[[2]]) %>% round(digits = 3)

cor(data$A[[2]], rr$ica$Mr[[1]]) %>% round(digits = 3)
cor(data$A[[2]], rr$ica$Mr[[2]]) %>% round(digits = 3)

modRV(data$A[[1]], rr$ica$Mr[[1]])
modRV(data$A[[2]], rr$ica$Mr[[2]])
modRV(data$A[[1]], rr$ica$Mr[[2]])
modRV(data$A[[2]], rr$ica$Mr[[1]])



cor(data$S[[1]], rr$ica$Sr[[1]]) %>% round(digits = 3)
cor(data$S[[2]], rr$ica$Sr[[2]]) %>% round(digits = 3)
cor(data$S[[3]], rr$ica$Sr[[3]]) %>% round(digits = 3)
modRV(data$S[[1]], rr$ica$Sr[[1]])
modRV(data$S[[2]], rr$ica$Sr[[2]])
modRV(data$S[[3]], rr$ica$Sr[[3]])

cor(data$A[[1]], rr$ica$Mr[[3]]) %>% round(digits = 3)
cor(data$A[[2]], rr$ica$Mr[[2]]) %>% round(digits = 3)
cor(data$A[[3]], rr$ica$Mr[[2]]) %>% round(digits = 3)
modRV(data$A[[1]], rr$ica$Mr[[3]])
modRV(data$A[[2]], rr$ica$Mr[[2]])
modRV(data$A[[3]], rr$ica$Mr[[2]])

