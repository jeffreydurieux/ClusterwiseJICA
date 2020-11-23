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

for(sim in 1:nrow(grid)){
  
  cat('\n')
  cat('\n')
  cat('------------')
  cat(sim)
  cat('------------')
  cat('\n')
  cat('\n')
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
  
  #### check tucker between signals of clusters####
  
  rhotuck <- numeric()
  combs <- combn(1:grid[sim,]$R, m = 2)
  for(com in 1:ncol(combs)){
    rhotuck[com] <- diag(Tucker(simdata$S[[combs[1,com]]],
                simdata$S[[combs[2,com]]])) %>% mean
  }
  combs <- rbind(combs,rhotuck)
  
  ##### analyse #######
  ptm <- proc.time()
  cjica <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
                           nc = grid[sim, ]$Q, starts = 100, scale = T)
  time <- proc.time() - ptm
  
  cjicaP <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
                           nc = grid[sim, ]$Q, starts = 1, scale = T,
                           rational = simdata$P )
  
  cjicaperbs <- list()
  for(perbs in 1:10){
    pert <- perturbation(p = simdata$P, percentage = 0.15)
    cjicaPert <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
                                 nc = grid[sim, ]$Q, starts = 1, scale = T,
                                 rational = pert)
    cjicaperbs[[perbs]] <- cjicaPert[[1]]
  }
  
  ##### evaluate ######
  loss100 <- sapply(seq_along(cjica), function(anom) tail(cjica[[anom]]$lossiter, n = 1))
  optimal <- cjica[[which.min(loss100)]]
  
  ### adjusted rand ###
  ari <- adjustedRandIndex(simdata$P, optimal$p)
  
  ### Tucker S ###
  tuckS <- list()
  for(tucks in 1:length(simdata$S) ){
    
    if(grid[sim,]$R == 2){
      tuckS[[tucks]] <- c(Tucker(simdata$S[[tucks]], optimal$ica$Sr[[1]]) %>% abs %>% 
        apply(MARGIN = 2, max) %>% mean, #i and 1
      Tucker(simdata$S[[tucks]], optimal$ica$Sr[[2]]) %>% abs %>% 
        apply(MARGIN = 2, max) %>% mean )#i and 1
    }else if(grid[sim,]$R == 3){
      tuckS[[tucks]] <- c(Tucker(simdata$S[[tucks]], optimal$ica$Sr[[1]]) %>% abs %>% 
        apply(MARGIN = 2, max) %>% mean, #i and 1
      Tucker(simdata$S[[tucks]], optimal$ica$Sr[[2]]) %>% abs %>% 
        apply(MARGIN = 2, max) %>% mean, #i and 2
      Tucker(simdata$S[[tucks]], optimal$ica$Sr[[3]]) %>% abs %>% 
        apply(MARGIN = 2, max) %>% mean )  #i and 3
      
    }else{
      tuckS[[tucks]] <- c(Tucker(simdata$S[[tucks]], optimal$ica$Sr[[1]]) %>% abs %>% 
        apply(MARGIN = 2, max) %>% mean , #i and 1
      Tucker(simdata$S[[tucks]], optimal$ica$Sr[[2]]) %>% abs %>% 
        apply(MARGIN = 2, max) %>% mean , #i and 2
      Tucker(simdata$S[[tucks]], optimal$ica$Sr[[3]]) %>% abs %>% 
        apply(MARGIN = 2, max) %>% mean , #i and 3
      Tucker(simdata$S[[tucks]], optimal$ica$Sr[[4]]) %>% abs %>% 
        apply(MARGIN = 2, max) %>% mean) #i and 4
    }
    
  }
  
  ### Tucker A ###
  #note: will cause error when no perfect P recovery. solve later after ALICE SIM
  # tuckA <- list()
  # for(tucks in 1:length(simdata$A) ){
  #   
  #   if(grid[sim,]$R == 2){
  #     tuckA[[tucks]] <- c(Tucker(simdata$A[[tucks]], optimal$ica$Mr[[1]]) %>% abs %>% 
  #                           apply(MARGIN = 2, max) %>% mean, #i and 1
  #                         Tucker(simdata$A[[tucks]], optimal$ica$Mr[[2]]) %>% abs %>% 
  #                           apply(MARGIN = 2, max) %>% mean )#i and 1
  #   }else if(grid[sim,]$R == 3){
  #     tuckA[[tucks]] <- c(Tucker(simdata$A[[tucks]], optimal$ica$Mr[[1]]) %>% abs %>% 
  #                           apply(MARGIN = 2, max) %>% mean, #i and 1
  #                         Tucker(simdata$A[[tucks]], optimal$ica$Mr[[2]]) %>% abs %>% 
  #                           apply(MARGIN = 2, max) %>% mean, #i and 2
  #                         Tucker(simdata$A[[tucks]], optimal$ica$Mr[[3]]) %>% abs %>% 
  #                           apply(MARGIN = 2, max) %>% mean )  #i and 3
  #     
  #   }else{
  #     tuckA[[tucks]] <- c(Tucker(simdata$A[[tucks]], optimal$ica$Mr[[1]]) %>% abs %>% 
  #                           apply(MARGIN = 2, max) %>% mean , #i and 1
  #                         Tucker(simdata$A[[tucks]], optimal$ica$Mr[[2]]) %>% abs %>% 
  #                           apply(MARGIN = 2, max) %>% mean , #i and 2
  #                         Tucker(simdata$A[[tucks]], optimal$ica$Mr[[3]]) %>% abs %>% 
  #                           apply(MARGIN = 2, max) %>% mean , #i and 3
  #                         Tucker(simdata$A[[tucks]], optimal$ica$Mr[[4]]) %>% abs %>% 
  #                           apply(MARGIN = 2, max) %>% mean) #i and 4
  #   }
  #   
  # }
  
  #### output list object #####
  output <- list()
  output$id <- grid[sim,]
  output$seed <- seed
  output$ari <- ari
  output$S <- tuckS
  output$Araw <- optimal$ica$Mr
  output$time <- time
  output$optimal <- optimal
  output$loss100 <- loss100
  output$rho <- combs
  output$cjicaP <- cjicaP[[1]]
  output$cjicaperbs <- cjicaperbs
  
  
  ext <- paste('/home/durieuxj/data/CJICA_results/Sim1/',
               'CJICA_sim1_',rows[sim], '.Rdata',sep = '')
  save(output,file = ext)
  
}


