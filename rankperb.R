# Tue Jun  8 14:28:57 2021
# Author: Jeffrey Durieux, MSc

rankperb <- function(Lir, nobj = 1){
  
  k <- sort(unique(Lir$newp))
  rank <- ssranking(Lir$ss)
  
  for(i in 1:nobj){
    kold <- Lir$newp[rank[i]]  
    ids <- which(kold != k)
    newm <- min(Lir$ss[rank[i], ids])
    knew <- which(Lir$ss[rank[i],] == newm)
    Lir$newp[rank[i]] <- knew
  }
  return(Lir$newp)
}