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

K <- 2; Q <- 2
data <- Simulate_CJICA(Nk = 20, Vm = 1000, K = K, Qm = Q, E = .01
                       , M = 2, type = 3)
str(data)

cjica <- ClusterwiseJICA(X = data$Xe, k = K, nc = Q, starts = 100)

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

Tucker(data$A[[1]], rr$ica$Mr[[2]]) %>% round(digits = 3)
Tucker(data$A[[2]], rr$ica$Mr[[1]]) %>% round(digits = 3)

modRV(data$A[[1]], rr$ica$Mr[[1]])
modRV(data$A[[2]], rr$ica$Mr[[2]])
modRV(data$A[[1]], rr$ica$Mr[[2]])
modRV(data$A[[2]], rr$ica$Mr[[1]])



Tucker(data$S[[1]], rr$ica$Sr[[3]]) %>% round(digits = 3)
Tucker(data$S[[2]], rr$ica$Sr[[2]]) %>% round(digits = 3)
Tucker(data$S[[3]], rr$ica$Sr[[1]]) %>% round(digits = 3)
modRV(data$S[[1]], rr$ica$Sr[[3]])
modRV(data$S[[2]], rr$ica$Sr[[2]])
modRV(data$S[[3]], rr$ica$Sr[[1]])

Tucker(data$A[[1]], rr$ica$Mr[[3]]) %>% round(digits = 3)
Tucker(data$A[[2]], rr$ica$Mr[[2]]) %>% round(digits = 3)
Tucker(data$A[[3]], rr$ica$Mr[[2]]) %>% round(digits = 3)
modRV(data$A[[1]], rr$ica$Mr[[3]])
modRV(data$A[[2]], rr$ica$Mr[[2]])
modRV(data$A[[3]], rr$ica$Mr[[2]])
