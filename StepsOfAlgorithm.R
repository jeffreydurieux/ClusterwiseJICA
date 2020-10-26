# Wed Oct 21 16:29:37 2020
# Author: Jeffrey Durieux, MSc


# What: steps of algorithm
# 1. generate P
# while no improvement in loss occur
# 2. extract ica on sorted data
# 3. compute Ahat, Xhat and Lir


source('dummydatascripts.R')
source('sortX.R')
source('ICAonList.R')
source('computeAhats.R')
source('computeXhats.R')
source('Avoid_nc_N.R')

ClusterwiseJICA <- function(X, k = 2, nc = 2, starts = 10){
  
  ResultsStarts <- list()
  
  for(start in 1:starts){
    
    totalSS <- sum(X^2)
    lossiter <- totalSS + 1
    iter <- 0
    
    repeat{
      iter <- iter + 1
      if(iter >= 2){
        cat('Iteration ', iter, ' loss value: ', lossiter[iter],'VAF:' ,Lir$vaf ,'\n')  
      }else{
        cat('Iteration ', iter, ' loss value: ', lossiter[iter],'\n') 
      }
      
      
      # algo step 1
      if(iter == 1){
        p <- CICA:::clusf(ncol(X), nClus = k)
      }else{
        p <- Lir$newp
      }
      List <- sortX(X, p)
      
      # algo step 2
      
      #### add stop warning over here ##### about nc <= k_n 
      icaparam <- ICAonList(List, nc = nc)
      
      # algo step 3
      Ahh <- Ahats(X = X, icapara = icaparam)
      Lir <- XhatsAndLir(X = X, Sr = icaparam$Sr, Ahats = Ahh)
      
      # avoid clus size lower than nc
      Lir$newp <- Avoid_nc_N(Lir$newp, Lir$lossvec, nc = nc)
      
      lossiter <- c(lossiter, Lir$loss)
      #abs(lossiter[iter + 1] - lossiter[iter])  < .00001
      if( lossiter[iter] - lossiter[iter + 1]  < .00001){
        break()
      }
    }
    
    
    out <- list()
    out$p <- Lir$newp
    out$ica <- icaparam
    out$lossiter <- lossiter
    ResultsStarts[[start]] <- out
    
    # make some sort of bubble sort algorithm so that only 2 starts are stored
    # output only the lowest start
    
  }
  
  
  
  return(ResultsStarts)
}
rm(res)

##### data contains 20 subjects with 2 underlying source components (two modalities each. since 2 underlying source components are used the number of clusters equals 2)

res <- ClusterwiseJICA(X = X, k = 4, nc = 3, starts = 10)


minloss <- sapply(seq_along(res), function(anom) tail(res[[anom]]$lossiter,1))
plot(minloss)

id <- minloss == min(minloss)

resmin <- res[id]

sols <- length(resmin)
id <- sample(sols,1)
rr <- resmin[[id]]
rr$p
k4 <- tail(rr$lossiter,n = 1)

table(rr$p)


library(CICA)
Tucker(SR1, rr$ica$Sr[[1]])
Tucker(SR2, rr$ica$Sr[[2]])

Tucker(SR1, rr$ica$Sr[[3]])
Tucker(SR2, rr$ica$Sr[[3]])


Tucker(A1, rr$ica$Mr[[2]])
Tucker(A2, rr$ica$Mr[[1]])
