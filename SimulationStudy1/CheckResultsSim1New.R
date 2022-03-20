# Tue Mar 15 17:45:45 2022
# Author: Jeffrey Durieux, MSc

# merge dfs from sim1 cjica

setwd('/Volumes/LaCie/backup_macbook/WeeklyBackup/P4_Data/SHARKDATA/Sim1DFS/')

library(gtools)

files <- dir()
files <- mixedsort(files)
files

load(files[1])

res <- df[,-c(1,2,3,4,5,6,7)]

for(i in 2:100){
  load(files[i])
  res <- res + df[,-c(1,2,3,4,5,6,7)]
}

results <- cbind(df[,c(1,2,3,4,5,6,7)], res)
range(results$aritrue)

splits <- split(1:3240, ceiling(seq_along(1:3240)/33))

load(files[101])

res2 <- df[,-c(1,2,3,4,5,6,7)]

id <- splits[[1]]
x <- get(load(files[100+1]))
res2 <- x[id,]
for(i in 2:99){
  id <- splits[[i]]
  x <- get(load(files[100+i]))
  res2 <- rbind(res2, x[id,])
  
}
results2 <- res2
#save(results2, file = '~/Downloads/sim1_N100_old.Rdata')

results3old <- rbind(results,results2old)
results3 <- rbind(results,results2)

save(results3, file = '~/Repo_temp/P4/SimulationStudy1/NewSim1df.Rdata')
load('~/Repo_temp/P4/SimulationStudy1/NewSim1df.Rdata')
library(doBy)

id <- results3$R==5 | results3$N ==50 | results3$Err==.6
sum(id)
res <- results3
res <- results3[!id,]

mean(results3$ari)
mean(res$ari)

summary_by(data = res, formula = ari~Q)
summary_by(data = res, formula = ari~R)
summary_by(data = res, formula = ari~N)
summary_by(data = res, formula = ari~rho)
summary_by(data = res, formula = ari~Err)
summary_by(data = res, formula = ari~lap) # what are the lap ra

mean(res$aritrue)
summary_by(data = res, formula = aritrue~Q)
summary_by(data = res, formula = aritrue~R)
summary_by(data = res, formula = aritrue~N)
summary_by(data = res, formula = aritrue~rho)
summary_by(data = res, formula = aritrue~lap) # what are the lap ranges


#m1
mean(res$arim1)
summary_by(data = res, formula = arim1~Q)
summary_by(data = res, formula = arim1~R)
summary_by(data = res, formula = arim1~N)
summary_by(data = res, formula = arim1~rho)
summary_by(data = res, formula = arim1~lap) # what are the lap ranges

#m2
mean(res$arim2)
summary_by(data = res, formula = arim2~Q)
summary_by(data = res, formula = arim2~R)
summary_by(data = res, formula = arim2~N)
summary_by(data = res, formula = arim2~rho)
summary_by(data = res, formula = arim2~lap) # what are the lap ra


ratio <- results3$N/results3$Q
mean(results3$ari)
mean(res$ari)

res <- cbind(res,ratio)
summary_by(data = res, formula = ari~ratio)


#### check loss

results3$freq_randomloss
tab <- table(results3$freq_randomloss)

id <- which(results3$freq_randomloss == 1)
rl1 <- results3[id,]

table(rl1$R)
rl11 <- rl1[which(rl1$R==3),]
mean(rl11$ari)
summary_by(data = rl11, formula = ari~Q)
summary_by(data = rl11, formula = ari~R)
summary_by(data = rl11, formula = ari~N)
summary_by(data = rl11, formula = ari~rho)
summary_by(data = rl11, formula = ari~lap)
rl11
rl11$randomloss - rl11$losstrue
rl11[68,]

### ALS probleem
bin <- ifelse(results3$freq_randomloss == 1, 1,0) # 0 is freq random loss > 1 . 1 is freq loss = 1
bindat <- cbind(results3,bin)
summary_by(data = bindat, ari~bin)
summary_by(data = bindat, ari~bin:Q)
summary_by(data = bindat, ari~bin:R)
summary_by(data = bindat, ari~bin:N)
summary_by(data = bindat, ari~bin:rho)
summary_by(data = bindat, ari~bin:lap)

#### smart start hcl
mean(results3$arikm)
summary_by(data = res, formula = arikm~Q)
summary_by(data = res, formula = arikm~R)
summary_by(data = res, formula = arikm~N)
summary_by(data = res, formula = arikm~rho)
summary_by(data = res, formula = arikm~lap) # what are the lap ra
