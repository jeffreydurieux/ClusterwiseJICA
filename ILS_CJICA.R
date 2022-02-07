# Tue Jun  8 11:21:40 2021
# Author: Jeffrey Durieux, MSc

# experimental ILS function for clusterwise JICA

# todo: right p is in the output but not the right Lir

ILSclusterwise <- function(X, p=NULL, Q, R, termination =20){
  
  start <- 1
  rational <- p
  k <- R
  nc <- Q
  
  # scaling, only for clusterwise joint ICA simulation design!
  # f1 <- sqrt(5000/sum(X[1:2500,]^2))
  # f2 <- sqrt(5000/sum(X[2501:5000,]^2))
  # X1 <- f1*X[1:2500,]
  # X2 <- f2*X[2501:5000,]
  # X <- rbind(X1,X2)
  
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
  
  #runs <- cbind(apply(lossvault,MARGIN = 2, adjustedRandIndex, TP),lossiter[-1])
  
  id <- which.min(lossiter[-1])
  
  out <- list()
  out$p <- lossvault[, id]
  out$loss_of_p <- min(lossiter)
  out$Lir <- Lir
  out$ICA <- icaparam
  out$lossiter <- lossiter
  
  return(out)
}
