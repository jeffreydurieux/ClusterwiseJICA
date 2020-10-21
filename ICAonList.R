# Wed Oct 21 09:46:37 2020
# Author: Jeffrey Durieux, MSc


# What: apply ICA on sorted data list with matrices
# algorithm step 2


ICAonList <- function(List, nc){
  Result <- lapply(List, FUN = ica::icafast, nc = nc) 
  S <- lapply(seq_along(Result), function(anom) Result[[anom]]$S)
  M <- lapply(seq_along(Result), function(anom) Result[[anom]]$M)
  
  params <- list(Sr=S,Mr=M)
  
  return(params)
}

#icaparam <- ICAonList(List, nc = 2)

