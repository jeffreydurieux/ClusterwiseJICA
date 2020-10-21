# Wed Oct 21 09:18:45 2020
# Author: Jeffrey Durieux, MSc

# dummy data for testing code quickly

# X => VxN 
# X = S %*% t(A)

library(ica)

set.seed(2407)

### create multimodal signals for cluster 1 and cluster 2

#### 2 multimodal signals R1
S1R1fun <- icasamp('b', nsamp = 1000, query = 'rnd')
S2R1fun <- icasamp('b', nsamp = 1000, query = 'rnd')

S1R1str <- icasamp('c', nsamp = 500, query = 'rnd')
S2R1str <- icasamp('c', nsamp = 500, query = 'rnd')

S1R1 <- c(S1R1fun,S1R1str)
S2R1 <- c(S2R1fun,S2R1str)

SR1 <- cbind(S1R1,S2R1)


#### 2 multimodal signals R2
S1R2fun <- icasamp('b', nsamp = 1000, query = 'rnd')
S2R2fun <- icasamp('b', nsamp = 1000, query = 'rnd')

S1R2str <- icasamp('c', nsamp = 500, query = 'rnd')
S2R2str <- icasamp('c', nsamp = 500, query = 'rnd')

S1R2 <- c(S1R2fun,S1R2str)
S2R2 <- c(S2R2fun,S2R2str)

SR2 <- cbind(S1R2,S2R2)

rm('S1R1fun', 'S1R1str','S1R2fun', 'S1R2str', 'S1R1', 'S2R1', 'S1R2', 'S2R2'
   ,'S2R1fun', 'S2R1str','S2R2fun', 'S2R2str')

n <- 20

A1 <- matrix(runif((n/2)*2, min = -2, max = 2), nrow = 10)
A2 <- matrix(runif((n/2)*2, min = -2, max = 2), nrow = 10)

A <- rbind(A1, A2)

X1 <- SR1 %*% t(A1)
X2 <- SR2 %*% t(A2)

X <- cbind(X1,X2)

### checks ICA ####
# Single ICA 
icawhole <- icafast(X, nc = 2)
CICA::Tucker(SR1 , icawhole$S)
CICA::Tucker(A, icawhole$M)

icar1 <- icafast(X[,1:10], nc = 2)
CICA::Tucker(SR1, icar1$S)
CICA::Tucker(A[1:10 ,], icar1$M)
#CICA::Tucker(A[1:10 ,], icar2$M)
CICA::Tucker(A1, icar1$M)

icar2 <- icafast(X[,11:20], nc = 2)
CICA::Tucker(SR2, icar2$S)
CICA::Tucker(SR2, icar1$S)
CICA::Tucker(A[11:20 ,], icar1$M)
CICA::Tucker(A[11:20 ,], icar2$M)
CICA::Tucker(A2, icar2$M)

rm('icar1', 'icar2', 'icawhole','S11','S12','S21','S22','S1','S2','Sa','Sb','n')



