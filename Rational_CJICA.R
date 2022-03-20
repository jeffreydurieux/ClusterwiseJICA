# Thu Mar 10 11:45:31 2022
# Author: Jeffrey Durieux, MSc

# What: rational start options for cjica
# try-out k-means, pam, hclust 
library(plotly)
library(cluster)

source('sortX.R')
source('ICAonList.R')
source('computeAhats.R')
source('computeXhats.R')
source('Avoid_nc_N.R')
source('Simulate_CJICA.R')
source('ClusterwiseJICA.R')
source('SearchEmptyClusters.R')
source('icafast_adjust.R')


Q <- c(2,5,10)
R <- c(2, 3, 4)
N <- c(20, 30, 50)
rho <- c(0, .50, .75)
E <- c(.2, .4, .75)
rep <- 1:20
grid <- expand.grid(Q = Q, R = R, N = N, rho = rho, E = E, rep = rep)



Nk = 50; Vm = 2500; K = 3; Qm = 2; E = .2; M =2; cor = .50; type = 4
Simulate_CJICA <- function(Nk, Vm, K=3, Qm, E, M){
  
  if(K != 3){
    stop('This simulation only applies to a K = 3 condition!')
  }
  
  # type 1: Sk %*% Ak
  # type 2: S %*% Ak
  # type 3: Sk %*% A
  # type 4: is type 1 but with pairwise correlated signals
  
  
  dnames <- c('b')
  
  P <- rep(1:K, each = Nk)
  
  
  if(type == 1 | type == 4){
    
    # for i in 1:2 --> modes
    # for k in 1:2 --> a vs bc and ab vs c
    
    Slist_m1 <- list()
    for(k in 1:2){
      
      S <- matrix(data = NA, nrow = Vm, ncol = Qm)
      for(q in 1:Qm){
        s <- icasamp(dname = sample(dnames, size = 1), query = 'rnd', nsamp = Vm)  
        S[, q] <- s
      }
      
      Slist_m1[[k]] <- S
    }
    
    Slist_m2 <- list()
    for(k in 1:2){
      
      S <- matrix(data = NA, nrow = Vm, ncol = Qm)
      for(q in 1:Qm){
        s <- icasamp(dname = sample(dnames, size = 1), query = 'rnd', nsamp = Vm)  
        S[, q] <- s
      }
      
      Slist_m2[[k]] <- S
    }
  }  
  
  
  
  Alist <- list()
  for(k in 1:3){
    
    A <- matrix(data = NA, nrow = Nk , ncol = Qm)
    for(q in 1:Qm){
      A[,q] <- runif(n = Nk, min = -2, max = 2)
    }
    Alist[[k]] <- A
  }
  
  
  S1 <- rbind(Slist_m1[[1]],Slist_m2[[1]])
  A <- S1 %*% t(Alist[[1]])
  
  S2 <- rbind(Slist_m1[[1]],Slist_m2[[2]])
  B <- S2 %*% t(Alist[[2]])
  
  S3 <- rbind(Slist_m1[[2]],Slist_m2[[2]])
  C <- S3 %*% t(Alist[[3]])
  
  X <- cbind(A,B,C)
  
  Xe <- addError(X, error = E)
  
  out <- list()
  out$Xe <- Xe
  out$X <- X
  out$P <- P
  out$S <- list(m1=Slist_m1, m2=Slist_m2)
  out$Scomb <- list(S1,S2,S3)
  out$A <- Alist
  return(out)  
}

set.seed(812)

ari_km <- numeric()
ari_pam <- numeric()
ari_hcl <- numeric()
for(i in 1:10){
  cat('loop number: ',i,'\n')
  
  test1 <- Simulate_CJICA(Nk = 20, Vm = 2500,K = 3, 
                          Qm = 2, E = .01, M = 2)
  
  X <- test1$Xe
  
  km_res <- kmeans(t(X), centers = 3, nstart = 100)
  ari_km[i] <- mclust::adjustedRandIndex(test1$P,km_res$cluster)
  
  pam_res <- pam(t(X), k = 3, diss = F)
  ari_pam[i] <- mclust::adjustedRandIndex(test1$P,pam_res$clustering)
  
  d <- dist(t(X))
  hcl <- hclust(d = d, method = 'ward.D2')
  plot(hcl)
  hcl_p <- cutree(hcl, k = 3)
  ari_hcl[i] <- mclust::adjustedRandIndex(test1$P, hcl_p)
}

mean(ari_km);mean(ari_pam);mean(ari_hcl)

#res <- ClusterwiseJICA(X = test1$Xe, k = 3, nc = 2, starts = 30, scale = T, verbose = T)

#losses <- sapply(seq_along(res), FUN = function(anom) tail(res[[anom]]$lossiter, n = 1) ) 
#plot(losses)
#ll <- sapply(seq_along(res), FUN = function(anom) tail(res[[anom]]$lossiter, n = 1) ) %>% which.min()
#res <- res[[ll]] 
#res$p
#mclust::adjustedRandIndex(test1$P,res$p)





library(multiway)
congru(test1$Scomb[[1]], res$ica$Sr[[1]])
congru(test1$Scomb[[2]], res$ica$Sr[[3]])
congru(test1$Scomb[[3]], res$ica$Sr[[2]])

#one modality
set.seed(812)
resm1 <- ClusterwiseJICA(X = test1$Xe[1:2500,], k = 2, nc = 2, starts = 30, scale = F, verbose = T)

losses <- sapply(seq_along(resm1), FUN = function(anom) tail(resm1[[anom]]$lossiter, n = 1) ) 
ll <- sapply(seq_along(resm1), FUN = function(anom) tail(resm1[[anom]]$lossiter, n = 1) ) %>% which.min()
resm1 <- resm1[[ll]] 
resm1$p
table(resm1$p)
mclust::adjustedRandIndex(test1$P,resm1$p)

#m2
set.seed(812)
resm2 <- ClusterwiseJICA(X = test1$Xe[2501:5000,], k = 3, nc = 2, starts = 30, scale = F, verbose = T)

losses <- sapply(seq_along(resm2), FUN = function(anom) tail(resm2[[anom]]$lossiter, n = 1) ) 
ll <- sapply(seq_along(resm2), FUN = function(anom) tail(resm2[[anom]]$lossiter, n = 1) ) %>% which.min()
resm2 <- resm2[[ll]] 
resm2$p
mclust::adjustedRandIndex(test1$P,resm2$p)