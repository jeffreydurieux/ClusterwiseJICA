# Thu Nov 19 14:22:52 2020
# Author: Jeffrey Durieux, MSc

# What: simulationscript for simulation 1 of clusterwise joint ica


#### Design ####

# 20 repetitions of each combination of the following factors:
#   - number of Q, 2, 5, 10
#   - number of R, 2, 3, 4
#   - N per cluster 20, 30, 50
#   - one-to-one correlation: no, low (.25), medium (.50)
#   - SNR  4, 1.5, .33 (20%, 40%, 75%)
# 
# Fixed: 5000 voxels, 100 random starts

# alice cluster info:
# will try to run it on the short partition.. Max time = 03:00:00
# maxjobsubmit is 100 so need to rewrite the code for --array=[1-100]


args <- commandArgs(TRUE)
#args <- as.numeric(args)

splits <- split(1:4860, ceiling(seq_along(1:4860)/49))
sp <- args[1]

rows <- splits[[sp]]


#### correlated clusters
setwd('/data/durieuxj/CJICACODE')
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



Q <- c(2,5,10)
R <- c(2, 3, 4)
N <- c(20, 30, 50)
rho <- c(0, .50, .75)
E <- c(.2, .4, .75)
rep <- 1:20
grid <- expand.grid(Q = Q, R = R, N = N, rho = rho, E = E, rep = rep)

grid <- grid[rows,]



