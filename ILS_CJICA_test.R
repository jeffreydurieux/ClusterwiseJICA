# Wed Mar  2 15:26:57 2022
# Author: Jeffrey Durieux, MSc
#https://link.springer.com/content/pdf/10.1007%2F978-3-319-07124-4_8.pdf

# ILS algorithm for CJICA

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

##### select data condition #####
sim <- which(grid$Q==2 & grid$R == 2 & grid$N == 50 & grid$rho ==0 & grid$E == 0.4)

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

rat <- c(rep(1,20), rep(2,20))#  , rep(3,20))
set.seed(1)
cjica <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
                         nc = grid[sim, ]$Q, starts = 10, scale = T)

losses <- sapply(seq_along(cjica), function(lam) tail(cjica[[lam]]$lossiter, n = 1) )
plot(losses)
opt <- which.min(losses)
min(losses)
cjica[[opt]]$p
mclust::adjustedRandIndex(simdata$P,cjica[[opt]]$p)
cjica[[opt]]$lossiter
x <- cjica[[opt]]


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


x <- cjica[[opt]]
X <- simdata$Xe;nc = grid[sim, ]$Q; k = grid[sim,]$R

f1 <- sqrt(5000/sum(X[1:2500,]^2))
f2 <- sqrt(5000/sum(X[2501:5000,]^2))
X1 <- f1*X[1:2500,]
X2 <- f2*X[2501:5000,]
X <- rbind(X1,X2)

X = simdata$Xe; k = grid[sim,]$R;
nc = grid[sim, ]$Q; starts = 1; scale = T
stepsize = 1

ILS_CJICA <- function(X, k, nc, scale = TRUE, iter, stepsize, titlevec){
  
  if(scale == TRUE){
    f1 <- sqrt(5000/sum(X[1:2500,]^2))
    f2 <- sqrt(5000/sum(X[2501:5000,]^2))
    X1 <- f1*X[1:2500,]
    X2 <- f2*X[2501:5000,]
    X <- rbind(X1,X2)
  }
  
  x <- ClusterwiseJICA(X = X, k = k, nc = nc, starts = 1, scale = scale)
  x <- x[[1]]
  
  lossvault <- x$Lir$loss
  pvault <- x$Lir$newp
  Lirvault <- list(x$Lir)
  Lir <- x$Lir
  
  n <- length(x$Lir$newp)
  it <- 0
  itvault <- 0
  losstrack <- x$Lir$loss
  temperature <- 1
  tempstep <- stepsize
  
  
  #it < iter & temperature < n
  while(it < iter){
    
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
      plot(losstrack, main = adjustedRandIndex(titlevec, newp))
      
      if(loss1 - loss2 < .01){
        
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

##### 
Tucker <- function(X, Y){
  return (diag(1 / sqrt(colSums(X^2))) %*% crossprod(X,Y) %*% diag(1 / sqrt(colSums(Y^2))) )
}
rm(test)
X <- simdata$Xe; k = length(simdata$S); nc = ncol(simdata$S[[1]])
test <- ILS_CJICA(X = X,k = k,nc = nc,scale = TRUE, iter = 1000,stepsize = 5)
test$lossvault
mclust::adjustedRandIndex(simdata$P, test$solution$p)

Tucker(simdata$S[[1]], test$solution$ica$Sr[[1]])
Tucker(simdata$S[[1]], test$solution$ica$Sr[[2]])

Tucker(simdata$S[[2]], test$solution$ica$Sr[[1]])
Tucker(simdata$S[[2]], test$solution$ica$Sr[[2]])



rat <- ClusterwiseJICA(X = X, k = 2, nc = nc, starts = 1, scale = TRUE, rational = simdata$P)
rat[[1]]$p
tail(rat[[1]]$lossiter,1)




