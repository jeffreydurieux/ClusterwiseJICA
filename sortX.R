# Wed Oct 21 09:18:01 2020
# Author: Jeffrey Durieux, MSc

# What: function that allocates subjects into matrix X based on partitioning
# Algorithm step 1

sortX <- function(X, p){
  ClusL <- length( unique(p) )
  NewList <- list()
  
  for(i in 1:ClusL){
    NewList[[i]] <- X[ ,p == i]
  }
  return(NewList)
}
p <- CICA:::clusf(20,2)
List <- sortX(X, p)


