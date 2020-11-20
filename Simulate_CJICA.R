# Mon Oct 26 13:42:45 2020
# Author: Jeffrey Durieux, MSc


# What: simulate clusterwise joint ICA data

# libraries used:
library(ica)

addError<-function(datablock,error)
{
  
  errorM<-replicate(ncol(datablock),rnorm(nrow(datablock)))
  errorM<-SSequal(errorM,datablock)
  errorlevel<-error/(1-error)
  
  res<-datablock + (errorM * sqrt(errorlevel))
  return(res)
}

SSequal<-function(m1,m2)
{
  res<-(m1/sqrt(sum(m1^2)) * sqrt(sum(m2^2))) #R
  return(res)
}


Simulate_CJICA <- function(Nk, Vm, K, Qm, E, M, type = 1, cor = .5){

  # type 1: Sk %*% Ak
  # type 2: S %*% Ak
  # type 3: Sk %*% A
  # type 4: is type 1 but with pairwise correlated signals
 
  
  dnames <- c('b')
  
  P <- rep(1:K, each = Nk)
  
  Slist <- list()
  if(type == 1 | type == 3 | type == 4){
    for(k in 1:K){
      
      S <- matrix(data = NA, nrow = Vm*M, ncol = Qm)
      for(q in 1:Qm){
        
        s <- numeric()
        for(m in 1:M){
          
          s <- c(s,s <- icasamp(dname = sample(dnames, size = 1),
                                query = 'rnd', nsamp = Vm)  )  
        }
        S[, q] <- s
      }
      Slist[[k]] <- S
      
      if(type == 4){
        
        
        if(K == 2){
          r <- matrix(c(1,cor,cor,1), nrow = 2)
          chol <- chol(r)  
          S2 <- matrix(data = NA, nrow = Vm*M, ncol = Qm)
          
          for(sig in 1:ncol(Slist[[1]])){
            ss <- cbind(Slist[[1]][,sig],
                        s <- icasamp(dname = sample(dnames, size = 1),
                                     query = 'rnd', nsamp = Vm*2)  )  
            ss <- ss %*% chol
            S2[,sig] <- ss[,2]
          }
          Slist[[2]] <- S2  
        }else if(K == 3){
          r <- matrix(c(1,cor,cor-.1,
                   cor,1,cor,
                   cor-.1, cor, 1),byrow = T, nrow = 3)
          
          
          chol <- chol(r)
          
          S2 <- matrix(data = NA, nrow = Vm*M, ncol = Qm)
          S3 <- matrix(data = NA, nrow = Vm*M, ncol = Qm)
          for(sig in 1:ncol(Slist[[1]])){
            ss <- cbind(Slist[[1]][,sig],
                        icasamp(dname = sample(dnames, size = 1),
                                     query = 'rnd', nsamp = Vm*2),
                        icasamp(dname = sample(dnames, size = 1),
                                query = 'rnd', nsamp = Vm*2))  
            ss <- ss %*% chol
            S2[,sig] <- ss[,2]
            S3[,sig] <- ss[,3]
          }
          Slist[[2]] <- S2
          Slist[[3]] <- S3
          
        }else if(K == 4){
          r <- matrix(c(1,cor,cor-.1, cor-.2,
                        cor,1, cor,cor-.1,
                        cor-.1,cor, 1,cor,
                        cor-.2,cor-.1, cor,1)
                      ,byrow = T, nrow = 4)
          
          chol <- chol(r)
          
          S2 <- matrix(data = NA, nrow = Vm*M, ncol = Qm)
          S3 <- matrix(data = NA, nrow = Vm*M, ncol = Qm)
          S4 <- matrix(data = NA, nrow = Vm*M, ncol = Qm)
          for(sig in 1:ncol(Slist[[1]])){
            ss <- cbind(Slist[[1]][,sig],
                        icasamp(dname = sample(dnames, size = 1),
                                     query = 'rnd', nsamp = Vm*2),
                        icasamp(dname = sample(dnames, size = 1),
                                query = 'rnd', nsamp = Vm*2),
                        icasamp(dname = sample(dnames, size = 1),
                                query = 'rnd', nsamp = Vm*2))  
            ss <- ss %*% chol
            S2[,sig] <- ss[,2]
            S3[,sig] <- ss[,3]
            S4[,sig] <- ss[,4]
          }
          Slist[[2]] <- S2
          Slist[[3]] <- S3
          Slist[[4]] <- S4
        }
        
      }
      
    }  
  }else if(type == 2){
    S <- matrix(data = NA, nrow = Vm*M, ncol = Qm)
    for(q in 1:Qm){
        s <- numeric()
        for(m in 1:M){
          s <- c(s,s <- icasamp(dname = sample(dnames, size = 1),
                                query = 'rnd', nsamp = Vm*2)  )  
          }
        S[, q] <- s
      }
      Slist[[1]] <- S
  }
  
  
  Alist <- list()
  for(k in 1:K){
    
    A <- matrix(data = NA, nrow = Nk, ncol = Qm)
    for(q in 1:Qm){
      A[,q] <- runif(n = Nk, min = -2, max = 2)
    }
    Alist[[k]] <- A
  }
  
  if(type == 1 | type == 4){
    X <- lapply(seq_along(Alist), function(anom)
      Slist[[anom]] %*% t(Alist[[anom]]))
  }else if(type == 2){
    X <- lapply(seq_along(Alist), function(anom)
      Slist[[1]] %*% t(Alist[[anom]]))
  }else{
    X <- lapply(seq_along(Alist), function(anom)
      Slist[[anom]] %*% t(Alist[[1]]))
  }
  
  X <- do.call(cbind, X)
  
  Xe <- addError(X, error = E)
  
  out <- list()
  out$Xe <- Xe
  out$X <- X
  out$P <- P
  out$S <- Slist
  out$A <- Alist
  return(out)  
}

# test1 <- Simulate_CJICA(Nk = 10, Vm = 2500,K = 4, 
#                         Qm = 5, E = .1, M = 2, type = 4, cor = .5)
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
