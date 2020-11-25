load('/Volumes/LaCie/MyData/CJICA/Sim2_1.Rdata')
RESULT

names <- rownames(RESULT)
ari <- paste(names, '_ARI', sep = '')
loss <- paste(names, '_LOSS', sep = '')
nm <- c(ari,loss)

tr <- t(RESULT)
res <- c(tr[1,],tr[2,])
names(res) <- nm
res

for(i in 2:243){
  ext <- paste('/Volumes/LaCie/MyData/CJICA/Sim2_',i,'.Rdata', sep = '')
  load(ext)
  tr <- t(RESULT)
  r <- c(tr[1,],tr[2,])
  res <- rbind(res, r)
}

Q <- c(2,5,10)
R <- c(2, 3, 4)
N <- c(20, 30, 50)
rho <- c(0, .50, .75)
E <- c(.2, .4, .75)
#rep <- 1:20
grid <- expand.grid(Q = Q, R = R, N = N, rho = rho, E = E)#, rep = rep)

res <- cbind(grid,res)
######## Quick statistics #######
library(doBy)

summary_by(data = res, formula = 
             cj_1_1_ARI+
             cj_1_2_ARI+
             cj_1_3_ARI+
             cj_1_4_ARI+
             cj_1_5_ARI ~
             rho*E*N, na.rm = T)

summary_by(data = res, formula = 
             cj_2_1_ARI+
             cj_2_2_ARI+
             cj_2_3_ARI+
             cj_2_4_ARI+
             cj_2_5_ARI ~
             ., na.rm = T)

summary_by(data = res, formula = 
             cj_3_1_ARI+
             cj_3_2_ARI+
             cj_3_3_ARI+
             cj_3_4_ARI+
             cj_3_5_ARI ~
             ., na.rm = T)


