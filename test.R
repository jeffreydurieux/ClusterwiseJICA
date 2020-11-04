# Mon Oct 26 14:28:14 2020
# Author: Jeffrey Durieux, MSc


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
load('~/Documents/datasets/Graz/AlzheimerfMRI_extra_aangevuld.RData')
X <- cbind(VBM4mm, FCdefaultModeNetwork4mm)
X <- t(X)
dim(X)
rm(correlations20components, correlations70components, partialCorrelations20components, partialCorrelations70components, FCdefaultModeNetwork2mm, FCdefaultModeNetwork4mm, FCexecutiveControl2mm, FCexecutiveControl4mm, VBM2mm, VBM4mm, ALFF2mm, ALFF4mm)

graztest <- ClusterwiseJICA(X = X, k = 2, nc = 5, scale = F, starts = 2)
graztest[[1]]$p

Alzheimer <- ifelse(Alzheimer == 0, 1, 2)
table(Alzheimer)

table(graztest[[1]]$p, Alzheimer)
graztest[[2]]$p
table(graztest[[2]]$p, Alzheimer)

tt <- table(graztest[[1]]$p, Alzheimer)
caret::confusionMatrix(tt)


K <- 2; Q <- 2; Nk = 50
data <- Simulate_CJICA(Nk = Nk, Vm = 2000, K = K, Qm = Q, E = .05
                       , M = 2, type = 1, cor = .71)
str(data)
cor(data$S[[1]],data$S[[2]])
cor(data$S[[1]],data$S[[3]])
#modRV(data$S[[1]],data$S[[2]])
time <- proc.time()
cjica <- ClusterwiseJICA(X = data$Xe, k = K, nc = Q, starts = 100)
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
Tucker(data$S[[1]], rr$ica$Sr[[1]]) %>% round(digits = 2)
Tucker(data$S[[2]], rr$ica$Sr[[2]]) %>% round(digits = 2)
Tucker(data$S[[1]], rr$ica$Sr[[2]]) %>% round(digits = 2)
Tucker(data$S[[2]], rr$ica$Sr[[1]]) %>% round(digits = 2)

modRV(data$S[[1]], rr$ica$Sr[[1]])
modRV(data$S[[2]], rr$ica$Sr[[2]])
modRV(data$S[[1]], rr$ica$Sr[[2]])
modRV(data$S[[2]], rr$ica$Sr[[1]])

Tucker(data$A[[1]], rr$ica$Mr[[1]]) %>% round(digits = 3)
Tucker(data$A[[1]], rr$ica$Mr[[2]]) %>% round(digits = 3)

Tucker(data$A[[2]], rr$ica$Mr[[1]]) %>% round(digits = 3)
Tucker(data$A[[2]], rr$ica$Mr[[2]]) %>% round(digits = 3)

modRV(data$A[[1]], rr$ica$Mr[[1]])
modRV(data$A[[2]], rr$ica$Mr[[2]])
modRV(data$A[[1]], rr$ica$Mr[[2]])
modRV(data$A[[2]], rr$ica$Mr[[1]])



Tucker(data$S[[1]], rr$ica$Sr[[1]]) %>% round(digits = 3)
Tucker(data$S[[2]], rr$ica$Sr[[2]]) %>% round(digits = 3)
Tucker(data$S[[3]], rr$ica$Sr[[3]]) %>% round(digits = 3)
modRV(data$S[[1]], rr$ica$Sr[[1]])
modRV(data$S[[2]], rr$ica$Sr[[2]])
modRV(data$S[[3]], rr$ica$Sr[[3]])

Tucker(data$A[[1]], rr$ica$Mr[[3]]) %>% round(digits = 3)
Tucker(data$A[[2]], rr$ica$Mr[[2]]) %>% round(digits = 3)
Tucker(data$A[[3]], rr$ica$Mr[[2]]) %>% round(digits = 3)
modRV(data$A[[1]], rr$ica$Mr[[3]])
modRV(data$A[[2]], rr$ica$Mr[[2]])
modRV(data$A[[3]], rr$ica$Mr[[2]])

