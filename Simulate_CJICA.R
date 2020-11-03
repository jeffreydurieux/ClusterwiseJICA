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
  # note type for only when K = 2
  
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
        r <- matrix(c(1,cor,cor,1), nrow = 2)
        chol <- chol(r)
        S2 <- matrix(data = NA, nrow = Vm*M, ncol = Qm)
        
        for(sig in 1:ncol(Slist[[1]])){
          ss <- cbind(Slist[[1]][,sig],
                      s <- icasamp(dname = sample(dnames, size = 1),
                                   query = 'rnd', nsamp = Vm)  )  
          ss <- ss %*% chol
          S2[,sig] <- ss[,2]
        }
        Slist[[2]] <- S2
      }
      
    }  
  }else if(type == 2){
    S <- matrix(data = NA, nrow = Vm*M, ncol = Qm)
    for(q in 1:Qm){
        s <- numeric()
        for(m in 1:M){
          s <- c(s,s <- icasamp(dname = sample(dnames, size = 1),
                                query = 'rnd', nsamp = Vm)  )  
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

test1 <- Simulate_CJICA(Nk = 10, Vm = 500,K = 2, 
                        Qm = 2, E = .1, M = 2, type = 1)
# str(test1)
# 
# test2 <- Simulate_CJICA(Nk = 10, Vm = 500,K = 2, 
#                         Qm = 2, E = .1, M = 2, type = 2)
# str(test2)
