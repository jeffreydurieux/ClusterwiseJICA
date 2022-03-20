# Wed Oct 21 09:46:37 2020
# Author: Jeffrey Durieux, MSc


# What: apply ICA on sorted data list with matrices
# algorithm step 2


ICAonList <- function(List, nc){
  
  #if else statement does not really matter. ica_adjust and ica::icafast give same results
  # included here only in testphase of coding
  #if(nc == 1){
    Result <- lapply(List, FUN = icafast_adjust, nc = nc) 
  #}else{
  #  Result <- lapply(List, FUN = ica::icafast, nc = nc) 
  #}

  S <- lapply(seq_along(Result), function(anom) Result[[anom]]$S)
  M <- lapply(seq_along(Result), function(anom) Result[[anom]]$M)
  
  params <- list(Sr=S,Mr=M)
  
  return(params)
}

#icaparam <- ICAonList(List, nc = 2)

