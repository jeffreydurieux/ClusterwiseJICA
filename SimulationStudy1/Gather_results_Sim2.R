# Wed Dec  2 13:25:49 2020
# Author: Jeffrey Durieux, MSc

# Gather results of clusterwise JICA simulation 2 (rational starts)

library(gtools)

setwd('/Volumes/LaCie/MyData/CJICA/Sim2/')
load('CJICA_sim2_1.Rdata')
RESULT

names <- rownames(RESULT)
ari <- paste(names, '_ARI', sep = '')
loss <- paste(names, '_LOSS', sep = '')
nm <- c(ari,loss)

tr <- t(RESULT)
res <- c(tr[1,],tr[2,])
names(res) <- nm
res

files <- dir()

files <- mixedsort(files)

files[1:10]

for(i in 2:4860){
  load(files[i])
  tr <- t(RESULT)
  r <- c(tr[1,],tr[2,])
  res <- rbind(res, r)
}

names <- c('f_SpecKernARI',
           'f_PamKernARI',
           'f_HclKernARI',
           'f_PamCosARI',
           'f_HclCosARI',
           
           'icaQ_SpecKernARI',
           'icaQ_PamKernARI',
           'icaQ_HclKernARI',
           'icaQ_PamCosARI',
           'icaQ_HclCosARI',
           
           'icaQR_SpecKernARI',
           'icaQR_PamKernARI',
           'icaQR_HclKernARI',
           'icaQR_PamCosARI',
           'icaQR_HclCosARI',
           
           'f_SpecKernLOSS',
           'f_PamKernLOSS',
           'f_HclKernLOSS',
           'f_PamCosLOSS',
           'f_HclCosLOSS',
           
           'icaQ_SpecKernLOSS',
           'icaQ_PamKernLOSS',
           'icaQ_HclKernLOSS',
           'icaQ_PamCosLOSS',
           'icaQ_HclCosLOSS',
           
           'icaQR_SpecKernLOSS',
           'icaQR_PamKernLOSS',
           'icaQR_HclKernLOSS',
           'icaQR_PamCosLOSS',
           'icaQR_HclCosLOSS'
           )

colnames(res) <- names

Q <- c(2,5,10)
R <- c(2, 3, 4)
N <- c(20, 30, 50)
rho <- c(0, .50, .75)
E <- c(.2, .4, .75)
rep <- 1:20
grid <- expand.grid(Q = Q, R = R, N = N, rho = rho, E = E, rep = rep)

dat <- cbind(grid,res)

save(dat, file = '/Volumes/LaCie/MyData/CJICA/Sim2/Sim2_results_df.Rdata')

