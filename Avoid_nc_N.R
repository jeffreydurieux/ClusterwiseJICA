# Thu Oct 22 11:02:22 2020
# Author: Jeffrey Durieux, MSc

# What: avoid nc < N

Avoid_nc_N <- function(newcluster, SSminVec, nc) {
  
  
  
  tab <- table(newcluster)
  
  if( all(tab >= nc) == F ){
    toofew <- which(tab < nc)
    
    # get number of elements need to pick per cluster
    topick <- nc - tab[toofew]
    topick <- as.integer(topick)
    
    id_to_pick_from <- which( !(newcluster %in% as.integer(toofew) ) )
    
    # sort based on loss value and get index
    ix <- sort.int(SSminVec, decreasing = T, index.return = T)
    id_to_add_to_small_clus <- ix$ix
    # remove ids small cluster
    id_to_add_to_small_clus <- id_to_add_to_small_clus[id_to_add_to_small_clus %in% id_to_pick_from]
    
    # loop over here
    #check here for multiple 
    for(pp in 1:length(topick)){
      worst_fit_pick <- id_to_add_to_small_clus[ 1:topick[pp] ]
      newcluster[worst_fit_pick] <- toofew[pp]
      id_to_add_to_small_clus <- id_to_add_to_small_clus[! id_to_add_to_small_clus %in% worst_fit_pick] 
    }
    
    tabb <- table(newcluster)
    if( any(tabb < nc) == TRUE){
      stop('Error: number of subject in a cluster is less than number of components. Even after attempts to reallocate worst fitting subjects. Select fewer components or fewer clusters')
    }
    
  }else{
    newcluster <- newcluster
  }
  
  
  return(newcluster)
}
# p <- c(1,1,1,1,1,1,1,1,
#        2,3,3,1,4,4,1,1,1,1,1,1)
# 
# p <- rep(1,240)
# p <- c(p, rep(2,10))
# table(p)
# SS <- runif(n = 250, min = .01, max = 5)
# Avoid_nc_N(newcluster = p, SSminVec = SS, nc = 20 ) %>% table()


# p <- c(1,1,1,1,1,1,1,1,
#        2,3,3,1,4,4,1,1,1,1,1,1)
# 
# p <- rep(1,240)
# p <- c(p, rep(1,10))
# table(p)
# SS <- runif(n = 250, min = .01, max = 5)
# Avoid_nc_N(newcluster = p, SSminVec = SS, nc = 20 ) %>% table()
