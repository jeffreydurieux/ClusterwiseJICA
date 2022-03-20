# Mon Mar 14 14:21:48 2022
# Author: Jeffrey Durieux, MSc

# what: new simulation with overlap and correlation structure cjica
# simulation run on Shark

#ext <- '/exports/fsw/durieuxj/P4/Sim1/'  

library(gtools)
library(mclust)
library(plotly)

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

#### Design ####

Q <- c(2,4,6)
R <- c(3, 5)
N <- c(50, 75)
rho <- c(0, .60, .94)
lap <- c(0,1,2)
Err <- c(.2, .4, .60)
rep <- 1:20

grid <- expand.grid(Q=Q, R=R, N=N, rho = rho, lap = lap, Err=Err, rep = rep)

grid <- grid[147,]

#args <- commandArgs(TRUE)
#args <- as.numeric(args)

#splits <- split(1:6480, ceiling(seq_along(1:6480)/65))
#sp <- args[1]

#rows <- splits[[sp]]

FindOptimalPermutSingle <- function( Sest , Strue, verbose = FALSE)
{
  # code to search the optimal permutation of estimated ICA components for
  # comparing it with simulated components
  # Author(s): Tom F. Wilderjans and minor adjustments by Jeffrey Durieux
  
  #JD: code from Tom, adjusted for matrix vs matrix comparison
  #Sest, Strue (nVoxels x nSources)
  
  library(gtools)
  N_sources = dim(Sest)[2]
  
  
  AllPerms = permutations( n = N_sources , r = N_sources , v = 1:N_sources )
  nPerms = dim(AllPerms)[1]
  
  #Find best permutation
  BestRecov = -9999
  BestPerm = -9999
  for( permtel in 1:nPerms )
  {
    if(verbose == TRUE)
    {
      if( (permtel%%50) == 0)
      {
        print( paste( "perm: " , permtel , "/" , nPerms ) )
      }
    }
    
    
    tp = AllPerms[permtel,]
    tempRecovBlock = matrix( -9999 , 1 , 1 )
    
    tempRecovBlock[1] = mean( abs( diag( Tucker(Strue ,
                                                Sest[, tp] ) ) ) )
    # niet nodig als het goed is
    tempRecov = mean(tempRecovBlock)
    
    if( permtel==1 )
    {
      BestRecov = tempRecov
      BestRecovBlock = tempRecovBlock
      BestPerm = tp
    }
    else
    {
      if( (tempRecov-BestRecov)>.0000000001 )
      {
        BestRecov = tempRecov
        BestRecovBlock = tempRecovBlock
        BestPerm = tp
      }
    }
    rm(tp,tempRecov,tempRecovBlock)
  }
  Out = list()
  Out$BestRecov = BestRecov
  Out$BestRecovBlock = BestRecovBlock
  Out$BestPerm = BestPerm
  Out$TuckerMatrix = Tucker(Strue , Sest[, BestPerm] )
  return(Out)
}

FindOptimalClusPermut <- function(Pest, Ptrue){
  # find optimal cluster permutation of estimated clustering
  # compared to simulated clustering
  # Author(s): Tom F. Wilderjans and minor adjustments by Jeffrey Durieux
  clus <- length(unique(Pest))
  
  AllPerms = gtools::permutations( n = clus , r = clus)
  nPerms = dim(AllPerms)[1]
  
  BestRecov = -9999
  BestPerm = -9999
  for( permtel in 1:nPerms )
  {
    if( (permtel%%50) == 0)
    {
      print( paste( "perm: " , permtel , "/" , nPerms ) )
    }
    
    tp = AllPerms[permtel,]
    tempRecovBlock = matrix( -9999 , 1 , 1 )
    
    tab <- table(Ptrue, Pest)
    
    tempRecovBlock[1] = sum( diag( tab[,tp] ) )
    
    tempRecov = mean(tempRecovBlock)
    
    if( permtel==1 )
    {
      BestRecov = tempRecov
      BestRecovBlock = tempRecovBlock
      BestPerm = tp
    }
    else
    {
      if( (tempRecov-BestRecov)>.0000000001 )
      {
        BestRecov = tempRecov
        BestRecovBlock = tempRecovBlock
        BestPerm = tp
      }
    }
    rm(tp,tempRecov,tempRecovBlock)
  }
  
  Out = list()
  Out$BestRecov = BestRecov
  Out$BestRecovBlock = BestRecovBlock
  Out$BestPerm = BestPerm
  return(Out)
  
}

Tucker <- function(X, Y){
  return (diag(1 / sqrt(colSums(X^2))) %*% crossprod(X,Y) %*% diag(1 / sqrt(colSums(Y^2))) )
}

modRV <- function(X, Y){
  
  if(nrow(X) != nrow(Y)){
    stop('Number of rows of input matrices are not equal')
  }
  
  XXtilde <- ( X %*% t(X) ) - diag (diag( X %*% t(X) ) )
  YYtilde <- ( Y %*% t(Y) ) - diag (diag( Y %*% t(Y) ) )
  
  res <-  ( t(c(XXtilde)) %*% c(YYtilde) ) /
    ( sqrt( ( t(c(XXtilde)) %*% c(XXtilde)) * ( t(c(YYtilde)) %*% c(YYtilde)) ) )
  
  
  return(res)
}

perturbation <- function(p, percentage = 0.1){
  
  clusters <- sort(unique(p))
  sel <- ceiling(length(p) * percentage )
  selected <- sample(1:length(p), size = sel, replace = F)
  
  if(length(selected) == 1){
    # change one cluster
    oriclus <- p[selected]
    newclus <- which(clusters != oriclus)
    
    if(length(newclus) > 1){
      newclus <- sample(newclus, size = 1)
    }
    
    np <- replace(p, selected, newclus)
    
  }else{
    # change multiple clusters
    np <- p
    for(i in 1:length(selected)){
      oriclus <- p[selected[i]]
      newclus <- which(clusters != oriclus)
      
      if(length(newclus) > 1){
        newclus <- sample(newclus, size = 1)
      }
      
      np <- replace(np, selected[i], newclus) # check if this works
    }
  }
  return(np)
}

clusf <- function(nBlocks, nClus) {
  #simplyfied cluster generation function using an equal probability
  clus <- GenerateRandomClustering(nBlocks, nClus, rep(c(1 / nClus), nClus))
  
  return(clus)
}


GenerateRandomClustering <- function(nElement , nClust , Prob = NULL)
{
  ####GenerateRandomClustering = for Random Starts
  
  # Author: Tom F. Wilderjans
  # nElement: number of elements to be clustered
  # nClust: number of clusters
  # Prob (1 x nClust): proportion of elements in each cluster
  
  # Added by Jeffrey Durieux: default Prob = equal cluster prob
  # This done to adjust code later on for potential cluster perbutation?
  
  if(is.null(Prob))
  {
    Prob <- rep(x = (1/nClust) , nClust)
  }
  
  
  BestClust = NULL
  ErrorEncountered = F
  
  if (!(length(Prob) == nClust))
  {
    cat('there should be as much probabilities as clusters')
    ErrorEncountered = T
  }
  
  if ((abs(sum(Prob) - 1) > .000000001) | (any(Prob < 0)))
  {
    cat('probabilities should sum to one (and cannot be negative)')
    ErrorEncountered = T
  }
  
  if (!(any(nClust == 1:nElement)))
  {
    cat("nClus should be a number between 1 and maximal number of datamatrices (length of DataList)")
    ErrorEncountered = T
  }
  
  if (!(ErrorEncountered))
  {
    if (nElement > nClust)
    {
      if (nClust == 1)
      {
        BestClust = rep(1 , times = nElement)
      }
      else
      {
        ProbVV = round(Prob * nElement)
        if (!(sum(ProbVV) == nElement) |
            (any(ProbVV < 1)))
          #not enough elements, or empty clusters
        {
          ProbVV = AdjustProb(ProbVV , nElement)
        }
        
        tempclus = rep(1:length(ProbVV) , ProbVV)
        BestClust = tempclus[sample(1:nElement,size = nElement,replace =
                                      FALSE)]
      }
    }
    else
    {
      BestClust = 1:nClust
    }
  }
  
  if (!(length(unique(BestClust)) == nClust))
  {
    BestClust = NULL
  }
  
  return(BestClust)
}


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




for(sim in 1:1){
  
  #cat('\n')
  #cat('\n')
  #cat('------------')
  #cat(sim)
  #cat('------------')
  #cat('\n')
  #cat('\n')
  seed <- as.numeric(rownames(grid))[sim]
  set.seed(seed)
  
  ##### simulate #####
  simdata <- Simulate_CJICA(Nk = grid[sim,]$N, 
                            Vm = 2500,
                            K = grid[sim, ]$R,
                            Qm = grid[sim, ]$Q,
                            Err = grid[sim, ]$Err,
                            M = 2,
                            cor = grid[sim, ]$rho,
                            con = grid[sim, ]$lap
  )
  
  ##### analyse #######
  #set.seed(parallel)
  ptm <- proc.time()
  cjica <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
                           nc = grid[sim, ]$Q, starts = 3000, scale = T)
  time <- proc.time() - ptm
  
  # cjicatrue <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
  #                              nc = grid[sim, ]$Q, starts = 1, scale = T,
  #                              rational = simdata$P )
  
  # cjicaperbs <- list()
  # for(perbs in 1:10){
  #   pert <- perturbation(p = simdata$P, percentage = 0.1)
  #   cjicaPert <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
  #                                nc = grid[sim, ]$Q, starts = 1, scale = T,
  #                                rational = pert)
  #   cjicaperbs[[perbs]] <- cjicaPert[[1]]
  # }
  # 
  # ##### cjica on mod one #####
  # f1 <- sqrt(5000/sum(simdata$Xe[1:2500,]^2))
  # X1 <- f1*simdata$Xe[1:2500,]
  # cjica_m1 <- ClusterwiseJICA(X = X1, k = grid[sim,]$R,
  #                             nc = grid[sim, ]$Q, starts = 100, scale = F)
  # 
  # ##### cjica on mod two #####
  # f2 <- sqrt(5000/sum(simdata$Xe[2501:5000,]^2))
  # X2 <- f2*simdata$Xe[2501:5000,]
  # cjica_m2 <- ClusterwiseJICA(X = X2, k = grid[sim,]$R,
  #                             nc = grid[sim, ]$Q, starts = 100, scale = F)
  # 
  # #### rational starts based on KM and HCL
  # km <- kmeans(x = t(simdata$Xe), centers = grid[sim, ]$R, nstart = 100)
  # 
  # cjica_km <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
  #                             nc = grid[sim, ]$Q, starts = 1, scale = T,
  #                             rational = km$cluster)
  # 
  # hcl <- hclust(dist( t(simdata$Xe)), method = 'ward.D2' )
  # cut <- cutree(hcl, k = grid[sim,]$R)
  # cjica_hcl <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
  #                              nc = grid[sim, ]$Q, starts = 1, scale = T,
  #                              rational = cut )
  # 
  # 
  ##### evaluate ######
  loss100 <- sapply(seq_along(cjica), function(anom) tail(cjica[[anom]]$lossiter, n = 1))
  optimal <- cjica[[which.min(loss100)]]
  
  # loss100m1 <- sapply(seq_along(cjica_m1), function(anom) tail(cjica_m1[[anom]]$lossiter, n = 1))
  # optimalm1 <- cjica_m1[[which.min(loss100m1)]]
  # 
  # loss100m2 <- sapply(seq_along(cjica_m2), function(anom) tail(cjica_m2[[anom]]$lossiter, n = 1))
  # optimalm2 <- cjica_m2[[which.min(loss100m2)]]
  # 
  ### adjusted rand ###
  ari <- adjustedRandIndex(simdata$P, optimal$p)
  # aritrueP <- adjustedRandIndex(simdata$P, cjicatrue[[1]]$p)
  # arim1 <- adjustedRandIndex(simdata$P, optimalm1$p)
  # arim2 <- adjustedRandIndex(simdata$P, optimalm2$p)
  
  
  ### Tucker S ###
  tucker_cor_lap <- unlist(TuckCheck(simdata$S))
  
  # add tucker congru 
  clusper <- FindOptimalClusPermut(optimal$p, simdata$P)
  
  if(grid[sim,]$R == 3){
    tucker_S <- mean(c(FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[1]]], Strue = simdata$S[[1]])$BestRecov,
                       FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[2]]], Strue = simdata$S[[2]])$BestRecov,
                       FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[3]]], Strue = simdata$S[[3]])$BestRecov))
    
  }else{
    tucker_S <- mean(c(FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[1]]], Strue = simdata$S[[1]])$BestRecov,
                       FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[2]]], Strue = simdata$S[[2]])$BestRecov,
                       FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[3]]], Strue = simdata$S[[3]])$BestRecov,
                       FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[4]]], Strue = simdata$S[[4]])$BestRecov,
                       FindOptimalPermutSingle(optimal$ica$Sr[[clusper$BestPerm[5]]], Strue = simdata$S[[5]])$BestRecov))
    
  }
  
  #### output list object #####
  output <- list()
  output$id <- grid[sim,]
  output$seed <- seed
  output$simdata <- simdata
  output$ari <- ari
  output$S <- tucker_S
  output$Araw <- optimal$ica$Mr
  output$time <- time
  output$optimal <- optimal
  output$loss100 <- loss100
  output$tuckercheck <- tucker_cor_lap
  #output$cjicatrue <- cjicatrue[[1]]
  #output$cjicaperbs <- cjicaperbs
  #output$m1 <- list(optimalm1, arim1)
  #output$m2 <- list(optimalm2, arim2)
  #output$km <- cjica_km
  #output$hcl <- cjica_hcl
  
  ext <- '~/Downloads/ExtraRandomStartsCJICA/'  
  ext <- paste(ext,'row_CJICA_sim1_',sim, '.Rdata',sep = '')
  save(output,file = ext)
  
}
load('~/Downloads/ExtraRandomStartsCJICA/row_254_CJICA_sim1_1.Rdata')
plot(output$loss100)
output$ari


#### ILS

###### ILS functions ######
scaleprob <- function(x){x/sum(x)}
uij2 <- function(x){sum(x^2)}
ssranking <- function(ss){
  ssscale <- apply(ss, 1, FUN = scaleprob)
  partcoef <- apply(ssscale, MARGIN = 2, FUN = uij2)
  sorted <- sort(partcoef, index.return=TRUE)
  return(sorted$ix)
}

rankperb <- function(Lir, nobj = 1){
  
  k <- sort(unique(Lir$newp))
  rank <- ssranking(Lir$ss)
  
  for(i in 1:nobj){
    kold <- Lir$newp[rank[i]]  
    ids <- which(kold != k)
    newm <- min(Lir$ss[rank[i], ids])
    knew <- which(Lir$ss[rank[i],] == newm)
    Lir$newp[rank[i]] <- knew
  }
  return(Lir$newp)
}





ILS_CJICA <- function(X, k, nc, scale = TRUE, iter, stepsize, titlevec, rational, maxtemp){
  
  if(scale == TRUE){
    f1 <- sqrt(5000/sum(X[1:2500,]^2))
    f2 <- sqrt(5000/sum(X[2501:5000,]^2))
    X1 <- f1*X[1:2500,]
    X2 <- f2*X[2501:5000,]
    X <- rbind(X1,X2)
  }
  
  x <- ClusterwiseJICA(X = X, k = k, nc = nc, starts = 1, scale = scale, rational = rational)
  x <- x[[1]]
  
  lossvault <- x$Lir$loss
  pvault <- x$Lir$newp
  Lirvault <- list(x$Lir)
  Lir <- x$Lir
  
  n <- maxtemp
  it <- 0
  itvault <- 0
  losstrack <- x$Lir$loss
  temperature <- 1
  tempstep <- stepsize
  
  
  #it < iter & temperature < n
  while(it < iter & temperature < n){
    
    it <- it + 1
    cat('Temperature equals :',temperature, fill = TRUE)
    cat('Iteration: ', it, fill = TRUE)
    newp <- rankperb(Lir = Lirvault[[which.min(lossvault)]], nobj = temperature)
    
    repeat{
      loss1 <- Lir$loss
      List <- sortX(X, newp)
      icaparam <- ICAonList(List, nc = nc)
      Ahh <- Ahats(X = X, icapara = icaparam)
      Lir <- XhatsAndLir(X = X, Sr = icaparam$Sr, Ahats = Ahh)
      loss2 <- Lir$loss
      newp <- Lir$newp
      
      loss1 - loss2
      losstrack <- c(losstrack,loss2)
      plot(losstrack, main = round(adjustedRandIndex(titlevec, newp), digits = 3))
      
      if(loss1 - loss2 < .00001){
        
        if(sign(loss1-loss2) == -1){
          #increase in loss
          cat('increase')
          temperature <- temperature + tempstep
          break()
        }else if(sign(loss1-loss2) == 0){
          # equal loss
          
          if(loss2 %in% lossvault){
            #if loss2 already in lossvault: increase temp
            temperature <- temperature + tempstep
            break()
          }
          
          if(loss2 > min(lossvault)){
            temperature <- temperature + tempstep
            break()
          }
          
          temperature <- 1
          lossvault <- c(lossvault, loss2)
          pvault <- cbind(pvault, newp)
          Lirvault <- c(Lirvault, list(Lir))
          itvault <- c(itvault,it)
          break()
        }
        
        break()
      } # end if
    }# end repeat
  }#end while
  
  out <- list()
  out$lossvault <- lossvault
  out$pvault <- pvault
  out$Lirvault <- Lirvault
  out$itvault <- itvault
  
  List <- sortX(X, newp)
  icaparam <- ICAonList(List, nc = nc)
  out$solution$ica <- icaparam
  out$solution$Lir <- Lir
  out$solution$p <- Lir$newp
  
  
  return(out)
}

set.seed(2407)
ILSres <- ILS_CJICA(X = simdata$Xe, k = grid$R, nc = grid$Q, scale = T, iter = 1000, stepsize = 1, 
                    titlevec = simdata$P, rational = output$optimal$p, maxtemp = 50)
min(ILSres$lossvault)
id <- which.min(ILSres$lossvault)

set.seed(2407)
test <- ClusterwiseJICA(X = simdata$Xe,k = grid$R, nc = grid$Q, scale = T,starts = 1, rational = perturbation(ILSres$pvault[,id], percentage = .1))
test[[1]]$lossiter


set.seed(2407)
ILSres2 <- ILS_CJICA(X = simdata$Xe, k = grid$R, nc = grid$Q, scale = T, iter = 100, stepsize = 1, 
                    titlevec = simdata$P, rational = test[[1]]$p, maxtemp = 100)
min(ILSres2$lossvault)
id <- which.min(ILSres2$lossvault)

set.seed(2407)
test2 <- ClusterwiseJICA(X = simdata$Xe,k = grid$R, nc = grid$Q, scale = T,starts = 1, rational = perturbation(ILSres2$pvault[,id], percentage = .1))
test2[[1]]$lossiter

set.seed(2407)
ILSres2 <- ILS_CJICA(X = simdata$Xe, k = grid$R, nc = grid$Q, scale = T, iter = 100, stepsize = 1, 
                     titlevec = simdata$P, rational = test2[[1]]$p, maxtemp = 15)
min(ILSres2$lossvault)
id <- which.min(ILSres2$lossvault)

test3 <- ClusterwiseJICA(X = simdata$Xe,k = grid$R, nc = grid$Q, scale = T,starts = 1, rational = perturbation(test3[[1]]$p, percentage = .1))
test3[[1]]$lossiter

