# Mon Mar 14 13:46:36 2022
# Author: Jeffrey Durieux, MSc

# Select values of overlap and correlation

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

#Nk = 50; Vm = 2500; K = 5; Qm = 2; Err = .2; M =2; cor = 0; type = 4; con = 2
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
    E <- S3 %*% t(Alist[[5]])
    
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

rho <- c(.94)
lap <- c(0, 1, 2)
clus <- c(3,5)
rep <- 1:10

pars <- expand.grid(rho=rho,lap = lap, clus=clus, rep)
cor <- numeric()
overlap <- numeric()
for(i in 1:nrow(pars)){
  sim <- Simulate_CJICA(Nk = 50, Vm = 2500, Qm = 2, Err = .1, M = 2, 
                        K = pars$clus[i], con = pars$lap[i], cor = pars$rho[i])  
  res <- TuckCheck(sim$S)
  cor[i] <- res[[1]]
  overlap[i] <- res[[2]]

}
df <- data.frame(pars, cor,overlap)

library(doBy)

corr <- summary_by(data = df, cor~rho:clus, FUN = c(mean,sd))
corr
lapp <- summary_by(data = df, overlap~lap:clus, FUN = c(mean,sd))
lapp

