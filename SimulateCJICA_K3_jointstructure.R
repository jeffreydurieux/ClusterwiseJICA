# Mon Mar  7 14:21:37 2022
# Author: Jeffrey Durieux, MSc


# What: simulation data gen test

setwd('~/Repo_temp/P4/EmpiricalAnalysis/')

source('../sortX.R')
source('../ICAonList.R')
source('../computeAhats.R')
source('../computeXhats.R')
source('../Avoid_nc_N.R')
source('../Simulate_CJICA.R')
source('../ClusterwiseJICA.R')
source('../SearchEmptyClusters.R')
source('../icafast_adjust.R')


# Q <- c(2,5,10)
# R <- c(2, 3, 4)
# N <- c(20, 30, 50)
# rho <- c(0, .50, .75)
# E <- c(.2, .4, .75)
# rep <- 1:20
# grid <- expand.grid(Q = Q, R = R, N = N, rho = rho, E = E, rep = rep)



Nk = 50; Vm = 2500; K = 5; Qm = 2; Err = .2; M =2; cor = 0; type = 4; con = 2
Simulate_CJICA <- function(Nk, Vm, K=3, Qm, Err, M, con, cor){
  
  if(K != 3  & K != 5){
    stop('This simulation only applies to a K = 3 condition!')
  }
  
  # type 1: Sk %*% Ak
  # type 2: S %*% Ak
  # type 3: Sk %*% A
  # type 4: is type 1 but with pairwise correlated signals
  
  
  dnames <- c('c')
  
  P <- rep(1:K, each = Nk)
  
     if(con != 0){
      Sbasem1 <- matrix(data = NA, nrow = Vm, ncol = Qm)
      for(q in 1:Qm){
        s <- icasamp(dname = sample(dnames, size = 1), query = 'rnd', nsamp = Vm)  
        Sbasem1[, q] <- s
      }
      
      
      r <- matrix(c(1,cor,
                    cor,1),byrow = T, nrow = 2)
      chol <- chol(r)
      
      Sbasem2 <- matrix(data = NA, nrow = Vm, ncol = Qm)
      
      for(sig in 1:ncol(Sbasem1)){
        ss <- cbind(Sbasem1[,sig],
                    icasamp(dname = sample(dnames, size = 1),
                            query = 'rnd', nsamp = Vm))  
        ss <- ss %*% chol
        Sbasem2[,sig] <- ss[,2]
      }
      
      #Tucker(Sbasem1,Sbasem2)
      
      Sbase <- rbind(Sbasem1,Sbasem2)
      
      Slist <- list()
      for(i in 1:K){
        Slist[[i]] <- replicate(n = Qm, 
                                runif(n = Vm*M, min = -con, max = con))
      }
      S <- lapply(seq_along(Slist), function(x) Sbase + Slist[[x]])
     }else{
       
       S <- list()
       for(clus in 1:K){
         Sbasem1 <- matrix(data = NA, nrow = Vm, ncol = Qm)
         for(q in 1:Qm){
           s <- icasamp(dname = sample(dnames, size = 1), query = 'rnd', nsamp = Vm)  
           Sbasem1[, q] <- s
         }
         
         
         r <- matrix(c(1,cor,
                       cor,1),byrow = T, nrow = 2)
         chol <- chol(r)
         
         Sbasem2 <- matrix(data = NA, nrow = Vm, ncol = Qm)
         
         for(sig in 1:ncol(Sbasem1)){
           ss <- cbind(Sbasem1[,sig],
                       icasamp(dname = sample(dnames, size = 1),
                               query = 'rnd', nsamp = Vm))  
           ss <- ss %*% chol
           Sbasem2[,sig] <- ss[,2]
         }
         
         
         S[[clus]] <- rbind(Sbasem1,Sbasem2)  
       }
       
    }
    
  
  
  Alist <- list()
  for(k in 1:K){
    
    A <- matrix(data = NA, nrow = Nk , ncol = Qm)
    for(q in 1:Qm){
      A[,q] <- runif(n = Nk, min = -2, max = 2)
    }
    Alist[[k]] <- A
  }
  
  if(K == 3){
    # A: 11   B:12    C:22
    
    S1 <- S[[1]]
    A <- S1 %*% t(Alist[[1]])
    
    S2 <- rbind(S[[1]][1:2500,],S[[2]][2501:5000,])
    B <- S2 %*% t(Alist[[2]])
    
    S3 <- S[[2]]
    C <- S3 %*% t(Alist[[3]])
    
    S <- list(S1,S2,S3)
    # Tucker(S1[1:2500,], S1[2501:5000,]) %>% round(digits = 4)
    # Tucker(S2[1:2500,], S2[2501:5000,]) %>% round(digits = 4)
    # Tucker(S3[1:2500,], S3[2501:5000,]) %>% round(digits = 4)
    # 
    # #combn(1:3,2)
    # Tucker(S1,S2)
    # Tucker(S1,S3)
    # Tucker(S2,S3)
    X <- cbind(A,B,C)
  }else if(K == 5){
    # A: 11   B:12    C:22      D:45  E:55
    S1 <- S[[1]]
    A <- S1 %*% t(Alist[[1]])
    
    S2 <- rbind(S[[1]][1:2500,],S[[2]][2501:5000,])
    B <- S2 %*% t(Alist[[2]])
    
    S3 <- S[[2]]
    C <- S3 %*% t(Alist[[3]])
    
    S4 <- S[[4]]
    D <- S4 %*% t(Alist[[4]])
    
    S5 <- rbind(S[[4]][1:2500,] , S[[5]][2501:5000,])
    E <- S5 %*% t(Alist[[5]]) #wasS3 first, a mistake
    
    S <- list(S1,S2,S3,S4,S5)
    
    # Tucker(S1[1:2500,], S1[2501:5000,]) %>% round(digits = 4)
    # Tucker(S2[1:2500,], S2[2501:5000,]) %>% round(digits = 4)
    # Tucker(S3[1:2500,], S3[2501:5000,]) %>% round(digits = 4)
    # Tucker(S4[1:2500,], S4[2501:5000,]) %>% round(digits = 4)
    # Tucker(S5[1:2500,], S5[2501:5000,]) %>% round(digits = 4)

    #combn(1:5,2)
    # Tucker(S1,S2)
    # Tucker(S1,S3)
    # Tucker(S1,S4)
    # Tucker(S1,S5)
    # 
    # Tucker(S2,S3)
    # Tucker(S2,S4)
    # Tucker(S2,S5)
    # 
    # Tucker(S3,S4)
    # Tucker(S3,S5)
    # 
    # Tucker(S4,S5)

    
    X <- cbind(A,B,C,D,E)
  }
  
  Xe <- addError(X, error = Err)
  
  out <- list()
  out$Xe <- Xe
  out$X <- X
  out$P <- P
  out$S <- S
  out$A <- Alist
  return(out)  
}

TuckCheck <- function(S){
  Tucker <- function(X, Y){
    return (diag(1 / sqrt(colSums(X^2))) %*% crossprod(X,Y) %*% diag(1 / sqrt(colSums(Y^2))) )
  }
  K <- length(S)
  
  if(K == 5){
    corMod <- mean(c(Tucker(S[[1]][1:2500,], S[[1]][2501:5000,]) %>% diag() %>% mean,
    Tucker(S[[2]][1:2500,], S[[2]][2501:5000,]) %>% diag() %>% mean,
    Tucker(S[[3]][1:2500,], S[[3]][2501:5000,]) %>% diag() %>% mean,
    Tucker(S[[4]][1:2500,], S[[4]][2501:5000,]) %>% diag() %>% mean,
    Tucker(S[[5]][1:2500,], S[[5]][2501:5000,]) %>% diag() %>% mean)
    )
    
    #combn(1:5,2)
    corClus <- mean(c(Tucker(S[[1]],S[[2]]) %>% diag() %>% mean,
    Tucker(S[[1]],S[[3]]) %>% diag() %>% mean,
    Tucker(S[[1]],S[[4]]) %>% diag() %>% mean,
    Tucker(S[[1]],S[[5]]) %>% diag() %>% mean,

    Tucker(S[[2]],S[[3]]) %>% diag() %>% mean,
    Tucker(S[[2]],S[[4]]) %>% diag() %>% mean,
    Tucker(S[[2]],S[[5]]) %>% diag() %>% mean,

    Tucker(S[[3]],S[[4]]) %>% diag() %>% mean,
    Tucker(S[[3]],S[[5]]) %>% diag() %>% mean,

    Tucker(S[[4]],S[[5]]) %>% diag() %>% mean
    ))
    
  }else if(K == 3){
    
    corMod <- mean(c(Tucker(S[[1]][1:2500,], S[[1]][2501:5000,]) %>% diag() %>% mean,
    Tucker(S[[2]][1:2500,], S[[2]][2501:5000,]) %>% diag() %>% mean,
    Tucker(S[[3]][1:2500,], S[[3]][2501:5000,]) %>% diag() %>% mean))

    #combn(1:3,2)
    corClus <- mean(c(Tucker(S[[1]],S[[2]])%>% diag() %>% mean,
    Tucker(S[[1]],S[[3]])%>% diag() %>% mean,
    Tucker(S[[2]],S[[3]])%>% diag() %>% mean))
  }
   
  return(list(corMod, corClus))
}

 # set.seed(812)
 # test1 <- Simulate_CJICA(Nk = 50, Vm = 2500, K = 5, 
 #                         Qm = 2, Err = .1, M = 2, con = 1.5, cor = 0)
 # TuckCheck(test1$S)
 # 
 # 
 # res <- ClusterwiseJICA(X = test1$Xe, k = 3, nc = 2, starts = 100, scale = T, verbose = T)
 # 
 # losses <- sapply(seq_along(res), FUN = function(anom) tail(res[[anom]]$lossiter, n = 1) ) 
 # plot(losses)
 # ll <- sapply(seq_along(res), FUN = function(anom) tail(res[[anom]]$lossiter, n = 1) ) %>% which.min()
 # res <- res[[ll]] 
 # res$p
 # mclust::adjustedRandIndex(test1$P,res$p)
 # 
 # library(multiway)
 # congru(test1$S[[1]], res$ica$Sr[[3]])
 # congru(test1$S[[2]], res$ica$Sr[[1]])
 # congru(test1$S[[3]], res$ica$Sr[[2]])
 # 
 # #one modality
 # set.seed(812)
 # resm1 <- ClusterwiseJICA(X = test1$Xe[1:2500,], k = 3, nc = 2, starts = 30, scale = F, verbose = T)
 # 
 # losses <- sapply(seq_along(resm1), FUN = function(anom) tail(resm1[[anom]]$lossiter, n = 1) ) 
 # plot(losses)
 # ll <- sapply(seq_along(resm1), FUN = function(anom) tail(resm1[[anom]]$lossiter, n = 1) ) %>% which.min()
 # resm1 <- resm1[[ll]] 
 # resm1$p
 # table(resm1$p)
 # mclust::adjustedRandIndex(test1$P,resm1$p)
 # 
 # #m2
 # set.seed(812)
 # resm2 <- ClusterwiseJICA(X = test1$Xe[2501:5000,], k = 3, nc = 2, starts = 30, scale = F, verbose = T)
 # 
 # losses <- sapply(seq_along(resm2), FUN = function(anom) tail(resm2[[anom]]$lossiter, n = 1) ) 
 # ll <- sapply(seq_along(resm2), FUN = function(anom) tail(resm2[[anom]]$lossiter, n = 1) ) %>% which.min()
 # resm2 <- resm2[[ll]] 
 # resm2$p
 # mclust::adjustedRandIndex(test1$P,resm2$p)
 #   
# str(test1)
# 
# combn(1:4,2)
# 
# mean(
# c(
# modRV(test1$S[[1]], test1$S[[2]]),
# modRV(test1$S[[1]], test1$S[[3]]),
# modRV(test1$S[[1]], test1$S[[4]]),
# modRV(test1$S[[2]], test1$S[[3]]),
# modRV(test1$S[[2]], test1$S[[4]]),
# modRV(test1$S[[3]], test1$S[[4]]) ) )
# 
# modRV(test1$S[[1]], test1$S[[2]])
# modRV(test1$S[[1]], test1$S[[3]])
# modRV(test1$S[[1]], test1$S[[4]])
# modRV(test1$S[[2]], test1$S[[3]])
# modRV(test1$S[[2]], test1$S[[4]])
# modRV(test1$S[[3]], test1$S[[4]]) 


# Tucker(test1$S[[1]], test1$S[[2]]) %>% round(digits = 3)
# Tucker(test1$S[[1]], test1$S[[3]]) %>% round(digits = 3)
# Tucker(test1$S[[1]], test1$S[[4]]) %>% round(digits = 3)
# Tucker(test1$S[[2]], test1$S[[3]]) %>% round(digits = 3)
# Tucker(test1$S[[2]], test1$S[[4]]) %>% round(digits = 3)
# Tucker(test1$S[[3]], test1$S[[4]]) %>% round(digits = 3)
# 
# 
# 
# test1 <- Simulate_CJICA(Nk = 10, Vm = 2500,K = 3,
#                        Qm = 2, E = .1, M = 2, type = 4)
#str(test1)

# combn(1:3,2)
# mean( c(
# modRV(test1$S[[1]], test1$S[[2]]) ,
# modRV(test1$S[[1]], test1$S[[3]]) ,
# modRV(test1$S[[2]], test1$S[[3]]) ))
# 
# modRV(test1$S[[1]], test1$S[[2]]) 
# modRV(test1$S[[1]], test1$S[[3]]) 
# modRV(test1$S[[2]], test1$S[[3]]) 
# 
# 
# test1 <- Simulate_CJICA(Nk = 10, Vm = 2500,K = 2, 
#                         Qm = 2, E = .1, M = 2, type = 4)
# str(test1)
# 
# combn(1:2,2)
# 
# modRV(test1$S[[1]], test1$S[[2]]) 
# 
# 
# 
# 
# Tucker(test1$S[[1]], test1$S[[2]])
# Tucker(test1$S[[2]], test1$S[[3]])
# Tucker(test1$S[[3]], test1$S[[4]])
# # 
# # test2 <- Simulate_CJICA(Nk = 10, Vm = 500,K = 2, 
# #                         Qm = 2, E = .1, M = 2, type = 2)
# # str(test2)

