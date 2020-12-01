SearchEmptyClusters <- function(nClus, newcluster, SSminVec) {
  
  OriCluster <- 1:nClus
  
  test <- sapply(OriCluster, FUN = '%in%', newcluster)
  
  #test result = no empty clusters so return original newcluster
  
  if ( all( test == TRUE) ){
    newcluster <- newcluster
  }else{
    
    EmptyClusters <- which(test == FALSE)
    singletonnames <- names(which( (table(newcluster)  == 1) == TRUE))
    singletons <- as.numeric(singletonnames)
    id <- which(newcluster %in% singletons == T)
    
    SSminVec[id] <- 0
    
    #worst <- sort( sapply( SSList, FUN = max), decreasing = TRUE)
    worst <- sort( SSminVec, decreasing = TRUE)
    
    #remove worst of singletons, otherwise empties will occur
    
    Index <- sapply( seq_along(EmptyClusters),
                     function(i) FUN = which( SSminVec == worst[i] ) )
    
    # if ties occur in SSminVec
    if( is.null(ncol(Index)) == FALSE ){
      Index <- Index[,1]
      Index <- sample(Index, size = length(EmptyClusters))
    }
    
    for(i in 1:length(Index)){
      newcluster <- replace(newcluster, Index[i], EmptyClusters[i])
      newcluster
    }
    
    
  }# else some emptyclusters
  if( length(unique(newcluster)) != nClus ){
    cat('SearchEmptyCluster, empty occurred')
  }
  return(newcluster)
}
#p <- rep(2,100)
#p <- c(p, rep(3,100))
#table(p)

#newp <- SearchEmptyClusters(nClus = 4, newcluster = p, SSminVec = Lir$lossvec)
#table(newp)
#Avoid_nc_N(newcluster = newp, SSminVec = Lir$lossvec, nc = 2) %>% table
