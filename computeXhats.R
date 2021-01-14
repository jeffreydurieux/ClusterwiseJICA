# Wed Oct 21 11:49:35 2020
# Author: Jeffrey Durieux, MSc


# What: computation of Xhats 
# algorithm step 3.2

# Xhat <- Sr %*% Air

XhatsAndLir <- function(X, Sr, Ahats){
  nClus <- length(Sr)
  N <- ncol(X)
  
  ss <- matrix(data = NA,nrow = N, ncol = nClus)
  for(outer in 1:N){
    
    for(inner in 1:nClus){
      xhat <- Sr[[inner]] %*% Ahats[[inner]][outer,]
      ss[outer,inner] <- sum( (X[,outer] - xhat)^2 )
    }
  }
  
  newp <- apply(ss, MARGIN = 1, which.min)
  lossvec <- apply(ss, MARGIN = 1, min)
  loss <- sum(lossvec)
  vaf <- ( sum(X^2)-loss) / sum(X^2)
  
  out <- list()
  out$newp <- newp
  out$lossvec <- lossvec
  out$loss <- loss
  out$vaf <- vaf
  out$ss <- ss
  return(out)
}

#Sr <- icaparam$Sr
#Ahats <- Ah

# pandloss <- XhatsAndLir(X,Sr,Ahats)
# pandloss$newp
# pandloss$lossvec
# pandloss$loss
# p
# newp
