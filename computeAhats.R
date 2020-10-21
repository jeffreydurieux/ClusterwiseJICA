# Wed Oct 21 11:38:53 2020
# Author: Jeffrey Durieux, MSc

# What: script to compute Ahats 
# algorithm step 3.1

Ahats <- function(X, icapara){

  nClus <- length(icapara$Sr)
  N <- ncol(X)
  
  AhatClus <- list()
  for(outer in 1:nClus){
    
    A <- matrix(data = NA, nrow = N, ncol(icapara$Mr[[1]]))
    
    for(inner in 1:N ){
      A[inner,] <- X[,inner] %*% icapara$Sr[[outer]] %*% 
        NMFN::mpinv( crossprod(icapara$Sr[[outer]]) )
    }
    AhatClus[[outer]] <- A
  }
  
  return(AhatClus)
}

Ah <- Ahats(X, icapara = icaparam)
