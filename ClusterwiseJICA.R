# Mon Oct 26 14:37:26 2020
# Author: Jeffrey Durieux, MSc

# Main function of clusterwise JICA
ClusterwiseJICA <- function(X, k = 2, nc = 2, starts = 10, scale = T){
  #### change this
  
  #if scale
  #f <- sqrt(1000/sum(X^2))
  #X <- f*X

  ResultsStarts <- list()
  
  for(start in 1:starts){
    
    totalSS <- sum(X^2)
    lossiter <- totalSS + 1
    iter <- 0
    
    repeat{
      iter <- iter + 1
      if(iter >= 2){
        cat('Start: ', start ,'Iteration ', iter, ' loss value: ', lossiter[iter],'VAF:' ,Lir$vaf ,'\n')  
      }else{
        cat('Start: ', start, 'Iteration ', iter, ' loss value: ', lossiter[iter],'\n') 
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
  }
  
  
  
  return(ResultsStarts)
}