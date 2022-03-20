# Tue Mar 15 09:46:40 2022
# Author: Jeffrey Durieux, MSc

# Script to extract results for new Simulation 1 for cjica 

library(gtools)
library(mclust)

args <- commandArgs(TRUE)
#args <- as.numeric(args)

splits <- split(1:3240, ceiling(seq_along(1:3240)/33))
sp <- args[1]

rows <- splits[[sp]]



#setwd('/exports/fsw/durieuxj/P4/Sim1/')
setwd('/Volumes/LaCie/backup_macbook/WeeklyBackup/P4_Data/SHARKDATA/Sim1/')

files <- dir()
id <- grep('N100', files)
files <- mixedsort(files[id])
files
Q <- c(2,4,6)
R <- c(3, 5)
N <- c(100) # c(50,75)
rho <- c(0, .60, .94)
lap <- c(0,1,2)
Err <- c(.2, .4, .60)
rep <- 1:20

grid <- expand.grid(Q=Q, R=R, N=N, rho = rho, lap = lap, Err=Err, rep = rep)


ari <- rep(0,3240)
randomloss <- rep(0,3240)
freq_randomloss <- rep(0,3240)
aritrue <- rep(0,3240)
losstrue <- rep(0,3240)
arim1 <- rep(0,3240)
arim2 <- rep(0,3240)
arikm <- rep(0,3240)
TuckerS <- rep(0,3240)
ariperb <- rep(0,3240)
lossperb <- rep(0,3240)
for(i in 1:length(grid)){
  load(files[rows[i]])
  
  # if(i == 1){
  #   id <- output$id  
  # }else{
  #   id <- rbind(id,output$id)
  # }
  # 
  
  id <- output$id  
  gridid <- as.numeric(rownames(id))
  ari[gridid] <- output$ari
  randomloss[gridid] <- min(output$loss100)
  freq_randomloss[gridid] <- length(which(output$loss100 == randomloss[gridid]))
  
  aritrue[gridid] <- adjustedRandIndex(output$cjicatrue$p, output$simdata$P)
  losstrue[gridid] <- tail(output$cjicatrue$lossiter, n = 1)
  
  arim1[gridid] <- adjustedRandIndex(output$m1[[1]]$p, output$simdata$P)
  arim2[gridid] <- adjustedRandIndex(output$m2[[1]]$p, output$simdata$P)
  
  arikm[gridid] <- adjustedRandIndex(output$hcl[[1]]$p, output$simdata$P)
  
  TuckerS[gridid] <- output$S
  idperb <- which.min(sapply(seq_along(output$cjicaperbs), function(anom) tail(output$cjicaperbs[[anom]]$lossiter, n = 1) ) )
  
  ariperb[gridid] <- adjustedRandIndex(output$cjicaperbs[[idperb]]$p, output$simdata$P)
  lossperb[gridid] <- min(sapply(seq_along(output$cjicaperbs), function(anom) tail(output$cjicaperbs[[anom]]$lossiter, n = 1) ) )
}

#df <- cbind(id,ari,randomloss,freq_randomloss, aritrue, losstrue, arim1, arim2, arikm, TuckerS, ariperb, lossperb)
#df2 <- cbind(ari,randomloss,freq_randomloss, aritrue, losstrue, arim1, arim2, arikm, TuckerS, ariperb, lossperb)
#s1 <- head(df1, n = 12)
#s2 <- head(df2, n = 12)

df <- cbind(grid,ari,randomloss,freq_randomloss, aritrue, losstrue, arim1, arim2, arikm, TuckerS, ariperb, lossperb)
ext <- paste('../Sim1ResultsDF_N100_',sp,'.Rdata', sep = '')
save(df, file = ext)