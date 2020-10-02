# Fri Oct  2 15:22:53 2020
# Author: Jeffrey Durieux, MSc

# What: script for simulating correlated signals

# Goal: simulate data according to clusterwise joint ica model (still need to formulate)
# this. Idea is to have....

# By applying a cholesky decomposition on a correlation matrix you can create a
# dataset with the specified correlation


#### example ####
# Cholesky decompsition --> rmat = R'R

r <- matrix(c(1,.3,.3,1), nrow = 2)
chol <- chol(r)
crossprod(chol)


So <- replicate(n = 2, expr = runif(n = 1000))
cor(So)

S <- So %*% chol
cor(S)


##### generate correlation matrix ######
r <- matrix(data = runif(3*3, min = .3, max = .5), nrow = 3)
r
#r[lower.tri(r)] <- 0
r[lower.tri(r)] <- t(r[upper.tri(r)])
diag(r) <- 1
r
isSymmetric(r)

S <- replicate(3,runif(1000, min = -2, max=2))
S <- S  %*% chol(r)
cor(S)

X <- S %*% matrix(data = rnorm(3*3), nrow = 3)

ica <- ica::icafast(X, nc = 3)
multiway::congru(ica$S, S)
