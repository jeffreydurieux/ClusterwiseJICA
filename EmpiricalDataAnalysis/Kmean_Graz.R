# Tue Jun  8 15:31:42 2021
# Author: Jeffrey Durieux, # Tue Jun  8 15:31:42 2021
# Author: Jeffrey Durieux, MSc

library(cluster)
library(factoextra)

# kmeans on graz data

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


km <- kmeans(x = X, centers = 2)
km$cluster

lab2 <- ifelse(lab==1,2,1)
tab <- table(rev(lab),km$cluster)
tab
caret::confusionMatrix(tab)

res <- clusGap(x = X, FUNcluster = kmeans,K.max = 5, B = 20, d.power = 1)
fviz_gap_stat(res)

avg_sil_values <- numeric()
for(i in 2:5){
  km.res <- kmeans(X, centers = i, nstart = 25)
  ss <- silhouette(km.res$cluster, dd)
  avg_sil_values[i] <- mean(ss[, 3])
  
}

k.values = 1:5

fviz_nbclust(X, FUNcluster = kmeans, method = 'silhouette')
