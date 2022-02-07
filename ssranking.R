# Tue Jun  8 14:30:01 2021
# Author: Jeffrey Durieux, MSc

ssranking <- function(ss){
  ssscale <- apply(ss, 1, FUN = scaleprob)
  #ssscale <- t(ssscale)
  partcoef <- apply(ssscale, MARGIN = 2, FUN = uij2)
  sorted <- sort(partcoef, index.return=TRUE)
  return(sorted$ix)
}