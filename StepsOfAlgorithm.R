# Wed Oct 21 16:29:37 2020
# Author: Jeffrey Durieux, MSc


# What: steps of algorithm
# 1. generate P
# while no improvement in loss occur
# 2. extract ica on sorted data
# 3. compute Ahat, Xhat and Lir


#### step 1 ####


ClusterwiseJICA <- function(X, k = 2, nc = 2, starts = 10){
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
    icaparam <- ICAonList(List, nc = nc)
    
    # algo step 3
    Ah <- Ahats(X = X, icapara = icaparam)
    Lir <- XhatsAndLir(X = X, Sr = icaparam$Sr, Ahats = Ah)
    
    lossiter <- c(lossiter, Lir$loss)
    #abs(lossiter[iter + 1] - lossiter[iter])  < .00001
    if( lossiter[iter] - lossiter[iter + 1]  < .00001){
      break()
    }
  }
  Lir$newp
  
  out <- list()
  out$p <- Lir$newp
  out$lossiter <- lossiter
  
  
  return(out)
}




