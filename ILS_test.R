# Thu Jan 14 12:12:09 2021
# Author: Jeffrey Durieux, MSc

# What: script for iterated local search for clusterwise procedures


library(mclust) #ARI
library(plotly) # '>' operator
source('sortX.R')
source('ICAonList.R')
source('computeAhats.R')
source('computeXhats.R')
source('Avoid_nc_N.R')
source('Simulate_CJICA.R')
source('ClusterwiseJICA.R')
source('SearchEmptyClusters.R')

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



Q <- c(2,5,10)
R <- c(2, 3, 4)
N <- c(20, 30, 50)
rho <- c(0, .50, .75)
E <- c(.2, .4, .75)

grid <- expand.grid(Q = Q, R = R, N = N, rho = rho, E = E)


sim <- which(grid$Q==5 & grid$R == 2 & grid$N == 20 & grid$rho ==0.5 & grid$E == 0.4)

# example Q=5 R =2 N = 20 rho = .5 E = 0.4 seed 110
seed <- as.numeric(rownames(grid))[sim]
set.seed(seed)

if(grid[sim, ]$rho == 0){
  type = 1
}else{
  type = 4
}

##### simulate #####
simdata <- Simulate_CJICA(Nk = grid[sim,]$N, 
                          Vm = 2500,
                          K = grid[sim, ]$R,
                          Qm = grid[sim, ]$Q,
                          E = grid[sim, ]$E,
                          M = 2,
                          cor = grid[sim, ]$rho,
                          type = type 
)

rat <- c(rep(1,20), rep(2,20), rep(3,20))
cjica <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
                         nc = grid[sim, ]$Q, starts = 10, scale = T)

losses <- sapply(seq_along(cjica), function(lam) tail(cjica[[lam]]$lossiter, n = 1) )
plot(losses)
opt <- which.min(losses)
min(losses)
cjica[[opt]]$p
mclust::adjustedRandIndex(simdata$P,cjica[[opt]]$p)


##### ssranking ####
scaleprob <- function(x){x/sum(x)}
uij2 <- function(x){sum(x^2)}
ssranking <- function(ss){
  ssscale <- apply(ss, 1, FUN = scaleprob)
  #ssscale <- t(ssscale)
  partcoef <- apply(ssscale, MARGIN = 2, FUN = uij2)
  sorted <- sort(partcoef, index.return=TRUE)
  return(sorted$ix)
}


##### rankperturbation ####
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


####### local opt with rank perp test #####

X <- simdata$Xe
scale=T
start <- 1
rational = cjica[[opt]]$p
#rational = lossvault[,id]
k = 2
nc = 5

if(scale == T){
  f1 <- sqrt(5000/sum(X[1:2500,]^2))
  f2 <- sqrt(5000/sum(X[2501:5000,]^2))
  X1 <- f1*X[1:2500,]
  X2 <- f2*X[2501:5000,]
  X <- rbind(X1,X2)
}




  
totalSS <- sum(X^2)
lossiter <- totalSS + 1
iter <- 0
  
rm(lossvault)
dynam <- 1
flagoccur <- 0
iter_since_lowest <- 0
lowestflag <- 0
##### manual loop #####

set.seed(seed)
while(iter_since_lowest <=20 & dynam <= 20 & lowestflag <= 20 ){
  iter <- iter + 1
  
  
  if(iter >= 2){
    cat('Start: ', start ,'Iteration ', iter, ' loss value: ', lossiter[iter],'VAF:' ,Lir$vaf ,'\n')  
  }else{
    cat('Start: ', start, 'Iteration ', iter, ' loss value: ', lossiter[iter],'\n') 
  }
  
  
  # algo step 1
  if(iter == 1){
    if(!is.null(rational)){
      p <- rational
      lossvault <- p
      t <- 0
      while( any( table(p)  < nc ) & t < 100 ){
        clusters <- 1:k
        
        id <- which(table(p) < nc)  
        id_to_take <- which(table(p) > nc)
        id_to_take <- which(p == id_to_take)
        
        s <- sample(id_to_take, size = 1)
        p[s] <- sample(id, size = 1)
        
        t <- t + 1
      }
      
      
    }else{
      p <- CICA:::clusf(ncol(X), nClus = k)
      lossvault <- p
    }
  }else{
    p <- Lir$newp
    lossvault <- cbind(lossvault,p)
  }
  List <- sortX(X, p)
  
  # algo step 2
  
  #### add stop warning over here ##### about nc <= k_n 
  icaparam <- ICAonList(List, nc = nc)
  
  # algo step 3
  Ahh <- Ahats(X = X, icapara = icaparam)
  Lir <- XhatsAndLir(X = X, Sr = icaparam$Sr, Ahats = Ahh)
  
  # avoid empty clusters
  if( length(unique(Lir$p)) < k ){
    Lir$newp <- SearchEmptyClusters(nClus = k, newcluster = Lir$newp, 
                                    SSminVec = Lir$lossvec)
  }
  
  # avoid clus size lower than nc
  Lir$newp <- Avoid_nc_N(Lir$newp, Lir$lossvec, nc = nc)
  
  lossiter <- c(lossiter, Lir$loss)
  
  print(lossiter)
  mclust::adjustedRandIndex(simdata$P,Lir$newp)
  if(length(lossiter) <= 20){
    title <- c('Random start ALS')  
  }else{
    title <- c('ILS based perturbation')
  }
  
  plot(lossiter[-1], ylim = c(2800,4000), ylab = 'Loss', main = title)
  #abs(lossiter[iter + 1] - lossiter[iter])  < .00001
  lossiter[iter] - lossiter[iter + 1]  < .00001
  lowest <- min(lossiter)
  lowestflag <- sum(lowest == lossiter)
  iter_since_lowest <- length(lossiter) - tail(which(lossiter==lowest), n = 1)
  if(lowestflag == 1 | lowestflag == 2 ){
    if(iter_since_lowest > 5){
      dynam <- dynam + 1
    }else{
      dynam <- 1
    }
  }else if(lowestflag == 4){
    if(iter_since_lowest > 5){
      dynam <- dynam + 1
    }else{
      dynam <- 2
    }
    flagoccur <- flagoccur + 1
  }else if(lowestflag == 6){
    if(iter_since_lowest > 5){
      dynam <- dynam + 1
    }else{
      dynam <- 3
    }
    flagoccur <- flagoccur + 1
  }else if(lowestflag == 8){
    if(iter_since_lowest > 5){
      dynam <- dynam + 1
    }else{
      dynam <- 4
    }
    flagoccur <- flagoccur + 1
  }
  # checks 
  Lir$newp
  
  #if climbing and equal to lowest occurred: increase rankperb
  # if not climbing but converged: increase rankperb by one
  # if only climbing: rankperb by 2
  # if bouncing above lowest (after 5 iters) increase dynam
  
  
  if(sign(lossiter[iter] - lossiter[iter+1]) == 0 & lossiter[iter+1] == lowest){
    dynam <- dynam + 1
    Lir$newp <- rankperb(Lir = Lir, nobj = dynam)
    cat('a \n')
  }else if(sign(lossiter[iter] - lossiter[iter+1]) == 0){
    Lir$newp <- rankperb(Lir = Lir, nobj = 1)  
    cat('b \n')
  }else if(sign(lossiter[iter] - lossiter[iter+1]) == -1){
    Lir$newp <- rankperb(Lir = Lir, nobj = 2)
    cat('c \n')
  }else if(lowestflag == 5){
    cat('Stop everything\n')
    print('stop everything')
  }
  
  if(iter_since_lowest > 5){
    dynam <- dynam + 1
    Lir$newp <- rankperb(Lir = Lir, nobj = dynam)
  }
  print(lossiter)
  print(mclust::adjustedRandIndex(simdata$P,Lir$newp))
  print(dynam)
  print(lowestflag)
  
}

  
  #check aris and loss 
  
  runs <- cbind(apply(lossvault,MARGIN = 2, adjustedRandIndex, simdata$P),lossiter[-1])
  id <- which.min(runs[,2])
  round(runs[id,],digits = 3)
  
  iter_since_lowest
  dynam
  lowestflag
  lossvault[,id]
  
  
  #### combine plots of random 10 start and ILS procedure
  # grid[110,] seed 110
  res <- c(losses[-1], lossiter[-1])
  plot(res)
  
  range(res)
  plot(losses[-1], xlim = c(0,128),ylim = c(2900,4000),col = 'darkgreen', 
       ylab = '', xlab = '')
  points(x = 10:128, y = res[-c(1:9)] )
  
  a1 <- locator(2)
  a2 <- locator(2)
  co.x <- cbind(a1$x,a2$x)
  co.y <- cbind(a1$y,a2$y)
  
  arrows(x0=co.x[2,], y0 = co.y[2,], x1 = co.x[1,], y1 = co.y[1,])
  text(x=co.x[2,], y=co.y[2,], labels=c("ARI = 0.14", "ARI = 1"), col=c("1", "1"), pos=c(1, 3), xpd=TRUE)
  mclust::adjustedRandIndex(simdata$P,cjica[[opt]]$p) # best of 10 random ARI
  round(runs[id,],digits = 3) # ILS ending ARI
  
  ##### extra info plot 
  a3 <- locator(2)
  a4 <- locator(2)
  co.x <- cbind(a3$x,a4$x)
  co.y <- cbind(a3$y,a4$y)
  
  arrows(x0=co.x[2,], y0 = co.y[2,], x1 = co.x[1,], y1 = co.y[1,])
  text(x=co.x[2,], y=co.y[2,], labels=c("Dynamische perturbatie", "Dynamische perturbatie"), col=c("1", "1"), pos=c(3, 3), xpd=TRUE)
  
  ####### notes ######
    
  # after convergence, add rankperb =1. 
  # if convergence equal stay rankperb ==1
  # if convergence 


#### function that does the ILS procedure #####
  
ILSclusterwise <- function(X, p=NULL, Q, R, termination =20){
  
  start <- 1
  rational <- p
  k <- R
  nc <- Q
  
  # scaling, only for clusterwise joint ICA simulation design!
  f1 <- sqrt(5000/sum(X[1:2500,]^2))
  f2 <- sqrt(5000/sum(X[2501:5000,]^2))
  X1 <- f1*X[1:2500,]
  X2 <- f2*X[2501:5000,]
  X <- rbind(X1,X2)
  
  totalSS <- sum(X^2)
  lossiter <- totalSS + 1
  iter <- 0
  
  dynam <- 1
  flagoccur <- 0
  iter_since_lowest <- 0
  lowestflag <- 0
  
  
  
  
  while(iter_since_lowest <=termination & dynam <= 20 & lowestflag <= termination){
    
    iter <- iter + 1
    
    
    if(iter >= 2){
      cat('Start: ', start ,'Iteration ', iter, ' loss value: ', lossiter[iter],'VAF:' ,Lir$vaf ,'\n')  
    }else{
      cat('Start: ', start, 'Iteration ', iter, ' loss value: ', lossiter[iter],'\n') 
    }
    
    
    # algo step 1
    if(iter == 1){
      if(!is.null(rational)){
        p <- rational
        lossvault <- p
        t <- 0
        while( any( table(p)  < nc ) & t < 100 ){
          clusters <- 1:k
          
          id <- which(table(p) < nc)  
          id_to_take <- which(table(p) > nc)
          id_to_take <- which(p == id_to_take)
          
          s <- sample(id_to_take, size = 1)
          p[s] <- sample(id, size = 1)
          
          t <- t + 1
        }
        
        
      }else{
        p <- CICA:::clusf(ncol(X), nClus = k)
        lossvault <- p
      }
    }else{
      p <- Lir$newp
      lossvault <- cbind(lossvault,p)
    }
    List <- sortX(X, p)
    
    # algo step 2
    
    #### add stop warning over here ##### about nc <= k_n 
    icaparam <- ICAonList(List, nc = nc)
    
    # algo step 3
    Ahh <- Ahats(X = X, icapara = icaparam)
    Lir <- XhatsAndLir(X = X, Sr = icaparam$Sr, Ahats = Ahh)
    
    # avoid empty clusters
    if( length(unique(Lir$p)) < k ){
      Lir$newp <- SearchEmptyClusters(nClus = k, newcluster = Lir$newp, 
                                      SSminVec = Lir$lossvec)
    }
    
    # avoid clus size lower than nc
    Lir$newp <- Avoid_nc_N(Lir$newp, Lir$lossvec, nc = nc)
    
    lossiter <- c(lossiter, Lir$loss)
    
    print(lossiter)
    #mclust::adjustedRandIndex(simdata$P,Lir$newp)
    plot(lossiter[-1])
    
    lossiter[iter] - lossiter[iter + 1]  < .00001
    
    lowest <- min(lossiter)
    lowestflag <- sum(lowest == lossiter)
    iter_since_lowest <- length(lossiter) - tail(which(lossiter==lowest), n = 1)
    
    if(lowestflag == 1 | lowestflag == 2 ){
      if(iter_since_lowest > 5){
        dynam <- dynam + 1
      }else{
        dynam <- 1
      }
    }else if(lowestflag == 4){
      if(iter_since_lowest > 5){
        dynam <- dynam + 1
      }else{
        dynam <- 2
      }
      flagoccur <- flagoccur + 1
    }else if(lowestflag == 6){
      if(iter_since_lowest > 5){
        dynam <- dynam + 1
      }else{
        dynam <- 3
      }
      flagoccur <- flagoccur + 1
    }else if(lowestflag == 8){
      if(iter_since_lowest > 5){
        dynam <- dynam + 1
      }else{
        dynam <- 4
      }
      flagoccur <- flagoccur + 1
    }
    
    
    #if climbing and equal to lowest occurred: increase rankperb
    # if not climbing but converged: increase rankperb by one
    # if only climbing: rankperb by 2
    # if bouncing above lowest (after 5 iters) increase dynam
    
    
    if(sign(lossiter[iter] - lossiter[iter+1]) == 0 & lossiter[iter+1] == lowest){
      dynam <- dynam + 1
      Lir$newp <- rankperb(Lir = Lir, nobj = dynam)
      #cat('a \n')
    }else if(sign(lossiter[iter] - lossiter[iter+1]) == 0){
      Lir$newp <- rankperb(Lir = Lir, nobj = 1)  
    }else if(sign(lossiter[iter] - lossiter[iter+1]) == -1){
      Lir$newp <- rankperb(Lir = Lir, nobj = 2)
    }
    
    if(iter_since_lowest > 5){
      dynam <- dynam + 1
      Lir$newp <- rankperb(Lir = Lir, nobj = dynam)
    }
    print(lossiter)
    

    
    
  }
  
  
  #check aris and loss 
  
  runs <- cbind(apply(lossvault,MARGIN = 2, adjustedRandIndex, simdata$P),lossiter[-1])
  id <- which.min(runs[,2])
  
  out <- list()
  out$p <- lossvault[, id]
  out$Lir <- Lir
  out$ICA <- icaparam
  
  return(out)
}

  
  ##### ILS function test #####
test <- ILSclusterwise(X = simdata$Xe, p = cjica[[opt]]$p, 
                       Q =grid[sim,]$Q ,R = grid[sim,]$R)
  
  
  
  
##### partitioning coefficient ######
test$p  
mclust::adjustedRandIndex(simdata$P,cjica[[opt]]$p)
mclust::adjustedRandIndex(simdata$P,test$p)
  
test$Lir$ss
SS <- cjica[[opt]]$Lir$ss

partitioningcoef <- function(SS){
  a <- apply(SS,MARGIN = 1, scaleprob)
  a <- 1-a
  
  partcoef <- apply(a, MARGIN = 2, FUN = max)
  return(sum(partcoef)/length(partcoef))
}
partitioningcoef(SS)

