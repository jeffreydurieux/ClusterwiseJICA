library(mclust) #ARI
library(plotly) # '>' operator
source('sortX.R')
source('ICAonList.R')
source('computeAhats.R')
source('computeXhats.R')
source('Avoid_nc_N.R')
source('Simulate_CJICA.R')
source('ClusterwiseJICA.R')
source('SearchEmptyClusters.R')

Q <- c(2,5,10)
R <- c(2, 3, 4)
N <- c(20, 30, 50)
rho <- c(0, .50, .75)
E <- c(.2, .4, .75)
rep <- 1:20
grid <- expand.grid(Q = Q, R = R, N = N, rho = rho, E = E, rep = rep)
sim = 809

cat('\n')
cat('\n')
cat('------------')
cat(sim)
cat('------------')
cat('\n')
cat('\n')
seed <- as.numeric(rownames(grid))[sim]
set.seed(seed)

if(grid[sim, ]$rho == 0){
  type = 1
}else{
  type = 4
}

##### simulate #####
simdata <- Simulate_CJICA(Nk = grid[sim,]$N, 
                          Vm = 2500,
                          K = grid[sim, ]$R,
                          Qm = grid[sim, ]$Q,
                          E = grid[sim, ]$E,
                          M = 2,
                          cor = grid[sim, ]$rho,
                          type = type 
)

#### check tucker between signals of clusters####


##### analyse #######

ptm <- proc.time()
cjica <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R,
                         nc = grid[sim, ]$Q, starts = 100, scale = T)
time <- proc.time() - ptm

