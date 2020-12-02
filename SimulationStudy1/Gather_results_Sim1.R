# Wed Nov 25 16:27:14 2020
# Author: Jeffrey Durieux, MSc

# Gather results CJICA simulation 1

library(gtools)
library(doBy)
library(stringr)
library(mclust)
source('~/Repositories/PhDthesis/P4/Simulate_CJICA.R')
setwd('/Volumes/LaCie/MyData/CJICA/Sim1/Sim1/')


Tucker <- function(X, Y){
  return (diag(1 / sqrt(colSums(X^2))) %*% crossprod(X,Y) %*% diag(1 / sqrt(colSums(Y^2))) )
}


files <- dir()
files <- mixedsort(files)
present <- files
present <- str_remove_all(present, "[CJICAsim.Rdata]")

present <- substring(present,4)
present
present <- as.numeric(present)
present

Q <- c(2,5,10)
R <- c(2, 3, 4)
N <- c(20, 30, 50)
rho <- c(0, .50, .75)
E <- c(.2, .4, .75)
rep <- 1:20
grid <- expand.grid(Q = Q, R = R, N = N, rho = rho, E = E, rep = rep)


grid <- grid[present,]

ARI <- numeric()
TIME <- numeric()
TuckS <- numeric()
LOSSrandom <- numeric()
LOSStrueP <- numeric()
LOSSperb <- numeric()
LOSSperbrange <- list()
LOSS100 <- list()
ARI_trueP <- numeric()
ARI_perbP <- numeric()

tmp <- proc.time()

# timing of this loop was: 1069 seconds/ 18 minutes
for(i in 1:nrow(grid)){
  cat('Extracting results from grid number: ',i,'\n')
  load(file = files[i])
  ARI[i] <- output$ari
  TIME[i] <- output$time[3]
  
  S <- output$S
  TuckS[i] <- sapply(S, FUN = max) %>% mean
  
  LOSSrandom[i] <- lossrandom <- tail(output$optimal$lossiter, n = 1)
  LOSStrueP[i] <- losstrueP <- tail(output$cjicaP$lossiter, n = 1)
  LOSSperb[i] <- lossperb <- min( sapply(seq_along(output$cjicaperbs), function(lam) 
    tail(output$cjicaperbs[[lam]]$lossiter, n = 1) ) )
  LOSSperbrange[[i]] <- range( sapply(seq_along(output$cjicaperbs), function(lam) 
    tail(output$cjicaperbs[[lam]]$lossiter, n = 1) ) )
  LOSS100[[i]] <- output$loss100
  
  seed <- as.numeric(rownames(grid))[i]
  set.seed(seed)
  
  if(grid[i, ]$rho == 0){
    type = 1
  }else{
    type = 4
  }
  
  ##### simulate #####
  simdata <- Simulate_CJICA(Nk = grid[i,]$N, 
                            Vm = 2500,
                            K = grid[i, ]$R,
                            Qm = grid[i, ]$Q,
                            E = grid[i, ]$E,
                            M = 2,
                            cor = grid[i, ]$rho,
                            type = type 
  )
  
  ARI_trueP[i] <- adjustedRandIndex(simdata$P, output$cjicaP$p)
  
  id_perb <- which.min (sapply(seq_along(output$cjicaperbs), 
         function(lam) tail(output$cjicaperbs[[lam]]$lossiter, n = 1)) )
  ARI_perbP[i] <- adjustedRandIndex(simdata$P,output$cjicaperbs[[id_perb]]$p)
  
  
}

gathering_time <- proc.time() - tmp

ARI <- numeric()
TIME <- numeric()
TuckS <- numeric()
LOSSrandom <- numeric()
LOSStrueP <- numeric()
LOSSperb <- numeric()
LOSSperbrange <- list()
LOSS100 <- list()
ARI_trueP <- numeric()
ARI_perbP <- numeric()

data <- cbind(grid, ARI, ARI_trueP, ARI_perbP, TuckS, LOSSrandom, LOSStrueP, LOSSperb)
save(data, file = '../Sim1_results_df.Rdata')
save(LOSS100, file = '../Sim1_results_LOSS100_list.Rdata')
save(LOSSperbrange, file = '../Sim1_results_LOSSperbrange_list.Rdata')


#### graph of loss 100 summary ####


