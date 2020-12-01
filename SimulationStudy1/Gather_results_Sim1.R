# Wed Nov 25 16:27:14 2020
# Author: Jeffrey Durieux, MSc

# Gather results CJICA simulation 1

library(gtools)
library(doBy)
library(stringr)

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
}



#### graph of loss 100 summary ####


library(HistogramTools)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

h1 <- range01(output1$loss100)
h2 <- range01(output2$loss100)
h1 <- hist(h1, plot = F)
h2 <- hist(h2, plot = F)
plot(h1)
plot(h2)
hadd <- AddHistograms(h1,h2)
par(mfrow=c(1,3))
plot(hadd)
# unlist matrix perbrange to matrix 
#matrix(unlist(LOSSperbrange, recursive = F), ncol = 2, byrow = T)