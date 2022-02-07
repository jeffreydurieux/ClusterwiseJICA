# Mon Oct 26 14:37:26 2020
# Author: Jeffrey Durieux, MSc

# Main function of clusterwise JICA
ClusterwiseJICA <- function(X, k = 2, nc = 2, starts = 10, scale = T, rational = NULL){
  #### change this
  
  if(scale == T){
    f1 <- sqrt(5000/sum(X[1:2500,]^2))
    f2 <- sqrt(5000/sum(X[2501:5000,]^2))
    X1 <- f1*X[1:2500,]
    X2 <- f2*X[2501:5000,]
    X <- rbind(X1,X2)
  }
  
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
        if(!is.null(rational)){
          p <- rational
          
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
          
        }
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
      
      # avoid empty clusters
      if( length(unique(Lir$p)) < k ){
        Lir$newp <- SearchEmptyClusters(nClus = k, newcluster = Lir$newp, 
                            SSminVec = Lir$lossvec)
      }
      
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
    out$Lir <- Lir
    ResultsStarts[[start]] <- out
  }
  
  
  
  return(ResultsStarts)
}
