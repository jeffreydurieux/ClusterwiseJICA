# Mon Nov 23 17:06:36 2020
# Author: Jeffrey Durieux, MSc


# smart start for clusterwise JICA
# use same grid as simulation study
# do pam, spectral clustering and hclust ward.D2 and complete

# do on full data, ICA data nc = Q per cluster or total Q (Q*R)

library(lsa) # cosine
library(cluster) # pam
library(mclust) # ari
library(kernlab) # spectral clustering
library(plotly) # '>' operator
source('sortX.R')
source('ICAonList.R')
source('computeAhats.R')
source('computeXhats.R')
source('Avoid_nc_N.R')
source('Simulate_CJICA.R')
source('ClusterwiseJICA.R')


Q <- c(2,5,10)
R <- c(2, 3, 4)
N <- c(20, 30, 50)
rho <- c(0, .50, .75)
E <- c(.2, .4, .75)
rep <- 1:20
grid <- expand.grid(Q = Q, R = R, N = N, rho = rho, E = E, rep = rep)


for(sim in 1:nrow(grid)){
  
  tmp <- proc.time()
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
  
  cos_sim_full <- cosine(simdata$Xe)
  ica1 <- icafast(simdata$Xe, nc = grid[sim,]$Q)
  ica2 <- icafast(simdata$Xe, nc = grid[sim,]$Q*grid[sim,]$R)
  Xhat1 <- ica1$S %*% t(ica1$M)
  Xhat2 <- ica2$S %*% t(ica2$M)
  cos_sim_ica1 <- cosine(Xhat1)
  cos_sim_ica2 <- cosine(Xhat2)
  
  ####### Part 1: full data ##########
  ### spectral clustering and kernelmatrix rbfdot 
  spec <- specc(x = t(simdata$Xe), centers = grid[sim,]$R, kernel = 'rbfdot')
  part1_p_spec_kern <- spec@.Data
  
  kernmat <- kernelMatrix(spec@kernelf@.Data, t(simdata$Xe))
  kernmat <- kernmat@.Data
  
  # pam and hcl with kernmat
  d <- as.dist(1 - kernmat)
  pam <- pam(x = d, k = grid[sim,]$R)
  part1_p_pam_kern <- pam$clustering
  
  hcl <- hclust(d = d, method = 'ward.D2')
  part1_p_hcl_kern <- cutree(hcl, k = grid[sim,]$R)
  
  
  # pam and hcl with cosine sim
  d <- as.dist(1 - cos_sim_full)
  pamcos <- pam(x = d, k = grid[sim,]$R)
  part1_p_pam_cos <- pamcos$clustering
  
  hcl <- hclust(d = d, method = 'ward.D2')
  part1_p_hcl_cos <- cutree(hcl, k = grid[sim,]$R)
  
  
  ####### Part 2: ica1 data ##########
  ### spectral clustering and kernelmatrix rbfdot 
  spec <- specc(x = t(Xhat1), centers = grid[sim,]$R, kernel = 'rbfdot')
  part2_p_spec_kern <- spec@.Data
  
  kernmat <- kernelMatrix(spec@kernelf@.Data, t(Xhat1))
  kernmat <- kernmat@.Data
  
  # pam and hcl with kernmat
  d <- as.dist(1 - kernmat)
  pam <- pam(x = d, k = grid[sim,]$R)
  part2_p_pam_kern <- pam$clustering
  
  hcl <- hclust(d = d, method = 'ward.D2')
  part2_p_hcl_kern <- cutree(hcl, k = grid[sim,]$R)
  
  # pam and hcl with cosine sim
  d <- as.dist(1 - cos_sim_ica1)
  pamcos <- pam(x = d, k = grid[sim,]$R)
  part2_p_pam_cos <- pamcos$clustering
  
  hcl <- hclust(d = d, method = 'ward.D2')
  part2_p_hcl_cos <- cutree(hcl, k = grid[sim,]$R)
  
  ####### Part 3: ica2 data ##########
  ### spectral clustering and kernelmatrix rbfdot 
  spec <- specc(x = t(Xhat2), centers = grid[sim,]$R, kernel = 'rbfdot')
  part3_p_spec_kern <- spec@.Data
  
  kernmat <- kernelMatrix(spec@kernelf@.Data, t(Xhat2))
  kernmat <- kernmat@.Data
  
  # pam and hcl with kernmat
  d <- as.dist(1 - kernmat)
  pam <- pam(x = d, k = grid[sim,]$R)
  part3_p_pam_kern <- pam$clustering
  
  hcl <- hclust(d = d, method = 'ward.D2')
  part3_p_hcl_kern <- cutree(hcl, k = grid[sim,]$R)
  
  # pam and hcl with cosine sim
  d <- as.dist(1 - cos_sim_ica2)
  pamcos <- pam(x = d, k = grid[sim,]$R)
  part3_p_pam_cos <- pamcos$clustering
  
  hcl <- hclust(d = d, method = 'ward.D2')
  part3_p_hcl_cos <- cutree(hcl, k = grid[sim,]$R)
  
  ######## part 4: CLJICA with smart starts ######
  
  
  cj_1_1 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part1_p_spec_kern, 
                            starts = 1)
  cj_1_2 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part1_p_pam_kern,
                            starts = 1)
  cj_1_3 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part1_p_hcl_kern,
                            starts = 1)
  cj_1_4 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part1_p_pam_cos,
                            starts = 1)
  cj_1_5 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part1_p_hcl_cos,
                            starts = 1)
  
  cj_2_1 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part2_p_spec_kern, 
                            starts = 1)
  cj_2_2 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part2_p_pam_kern,
                            starts = 1)
  cj_2_3 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part2_p_hcl_kern,
                            starts = 1)
  cj_2_4 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part2_p_pam_cos,
                            starts = 1)
  cj_2_5 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part2_p_hcl_cos,
                            starts = 1)
  
  cj_3_1 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part3_p_spec_kern, 
                            starts = 1)
  cj_3_2 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part3_p_pam_kern,
                            starts = 1)
  cj_3_3 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part3_p_hcl_kern,
                            starts = 1)
  cj_3_4 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part3_p_pam_cos,
                            starts = 1)
  cj_3_5 <- ClusterwiseJICA(X = simdata$Xe, k = grid[sim,]$R, 
                            nc = grid[sim,]$Q, rational = part3_p_hcl_cos,
                            starts = 1)
  
  ####### results ########
  # ari and loss
  ARI <- c(adjustedRandIndex(simdata$P, cj_1_1[[1]]$p),
           adjustedRandIndex(simdata$P, cj_1_2[[1]]$p),
           adjustedRandIndex(simdata$P, cj_1_3[[1]]$p),
           adjustedRandIndex(simdata$P, cj_1_4[[1]]$p),
           adjustedRandIndex(simdata$P, cj_1_5[[1]]$p),
           
           adjustedRandIndex(simdata$P, cj_2_1[[1]]$p),
           adjustedRandIndex(simdata$P, cj_2_2[[1]]$p),
           adjustedRandIndex(simdata$P, cj_2_3[[1]]$p),
           adjustedRandIndex(simdata$P, cj_2_4[[1]]$p),
           adjustedRandIndex(simdata$P, cj_2_5[[1]]$p),
           
           adjustedRandIndex(simdata$P, cj_3_1[[1]]$p),
           adjustedRandIndex(simdata$P, cj_3_2[[1]]$p),
           adjustedRandIndex(simdata$P, cj_3_3[[1]]$p),
           adjustedRandIndex(simdata$P, cj_3_4[[1]]$p),
           adjustedRandIndex(simdata$P, cj_3_5[[1]]$p))
  
  LOSS <- c(tail( cj_1_1[[1]]$lossiter, n = 1),
            tail( cj_1_2[[1]]$lossiter, n = 1),
            tail( cj_1_3[[1]]$lossiter, n = 1),
            tail( cj_1_4[[1]]$lossiter, n = 1),
            tail( cj_1_5[[1]]$lossiter, n = 1),
            
            tail( cj_2_1[[1]]$lossiter, n = 1),
            tail( cj_2_2[[1]]$lossiter, n = 1),
            tail( cj_2_3[[1]]$lossiter, n = 1),
            tail( cj_2_4[[1]]$lossiter, n = 1),
            tail( cj_2_5[[1]]$lossiter, n = 1),
            
            tail( cj_3_1[[1]]$lossiter, n = 1),
            tail( cj_3_2[[1]]$lossiter, n = 1),
            tail( cj_3_3[[1]]$lossiter, n = 1),
            tail( cj_3_4[[1]]$lossiter, n = 1),
            tail( cj_3_5[[1]]$lossiter, n = 1))
  
  names <- c('cj_1_1',
             'cj_1_2',
             'cj_1_3',
             'cj_1_4',
             'cj_1_5',
             
             'cj_2_1',
             'cj_2_2',
             'cj_2_3',
             'cj_2_4',
             'cj_2_5',
             
             'cj_3_1',
             'cj_3_2',
             'cj_3_3',
             'cj_3_4',
             'cj_3_5')
  
  RESULT <- data.frame(ARI,LOSS)
  rownames(RESULT) <- names
  time <- proc.time() - tmp
  time
  ext <- paste('Sim2_',grid[sim, ]$rep,'.Rdata',sep = '')
  
  save(RESULT, file = ext)
}