# Wed Dec  2 16:35:12 2020
# Author: Jeffrey Durieux, MSc

# What: Simulation 1 histograms of loss function values

load('/Volumes/LaCie/MyData/CJICA/Sim1/Sim1_results_LOSS100_list.Rdata')

library(HistogramTools)


Q <- c(2,5,10)
R <- c(2, 3, 4)
N <- c(20, 30, 50)
rho <- c(0, .50, .75)
E <- c(.2, .4, .75)
rep <- 1:20
grid <- expand.grid(Q = Q, R = R, N = N, rho = rho, E = E, rep = rep)


### example first dataset loss100 versus last
par(mfrow=c(1,2))
hist(LOSS100[[1]])
hist(LOSS100[[4860]])

library(HistogramTools)

range01 <- function(x){(x-min(x))/(max(x)-min(x))}

AddHisto <- function(xlist){
  len <- length(xlist)
  x1 <- range01(xlist[[1]])
  x1 <- hist(x1, plot = F)
  x2 <- range01(xlist[[2]])
  x2 <- hist(x2, plot = F)
  
  hadd <- AddHistograms(x1,x2)
  for(i in 3:len){
    x2 <- range01(xlist[[i]]) 
    x2 <- hist(x2, plot = F)
    hadd <- AddHistograms(hadd,x2)
  }
  return(hadd)
}

par(mfrow=c(1,1))
overall <- AddHisto(xlist = LOSS100)
plot(overall, main ='', xlab ='Overall summed histograms')


##### Q #####
par(mfrow=c(1,3))
id <- which(grid$Q == 2)
Q2 <- AddHisto(xlist = LOSS100[id])
plot(Q2, main ='', xlab ='Q2 summed histograms', ylim = c(0,50000))

id <- which(grid$Q == 5)
Q5 <- AddHisto(xlist = LOSS100[id])
plot(Q5, main ='', xlab ='Q5 summed histograms', ylim = c(0,50000))

id <- which(grid$Q == 10)
Q10 <- AddHisto(xlist = LOSS100[id])
plot(Q10, main ='', xlab ='Q10 summed histograms', ylim = c(0,50000))


##### R #####
par(mfrow=c(1,3))
id <- which(grid$R == 2)
R2 <- AddHisto(xlist = LOSS100[id])
plot(R2, main ='', xlab ='R2 summed histograms', ylim = c(0,40000))

id <- which(grid$R == 3)
R3 <- AddHisto(xlist = LOSS100[id])
plot(R3, main ='', xlab ='R3 summed histograms', ylim = c(0,40000))

id <- which(grid$R == 4)
R4 <- AddHisto(xlist = LOSS100[id])
plot(R4, main ='', xlab ='R4 summed histograms', ylim = c(0,40000))


##### N #####
par(mfrow=c(1,3))
id <- which(grid$N == 20)
N20 <- AddHisto(xlist = LOSS100[id])
plot(N20, main ='', xlab ='N20 summed histograms', ylim = c(0,30000))

id <- which(grid$N == 30)
N30 <- AddHisto(xlist = LOSS100[id])
plot(R3, main ='', xlab ='N30 summed histograms', ylim = c(0,30000))

id <- which(grid$N == 50)
N50 <- AddHisto(xlist = LOSS100[id])
plot(N50, main ='', xlab ='N50 summed histograms', ylim = c(0,40000))


##### rho #####
par(mfrow=c(1,3))
id <- which(grid$rho == 0)
rho0 <- AddHisto(xlist = LOSS100[id])
plot(rho0, main ='', xlab ='rho0 summed histograms', ylim = c(0,30000))

id <- which(grid$rho == 0.5)
rho0.5 <- AddHisto(xlist = LOSS100[id])
plot(rho0.5, main ='', xlab ='rho0.5 summed histograms', ylim = c(0,30000))

id <- which(grid$rho == 0.75)
rho0.75 <- AddHisto(xlist = LOSS100[id])
plot(rho0.75, main ='', xlab ='rho0.75 summed histograms', ylim = c(0,40000))


##### E #####
par(mfrow=c(1,3))
id <- which(grid$E == 0.20)
E20 <- AddHisto(xlist = LOSS100[id])
plot(E20, main ='', xlab ='E20 summed histograms', ylim = c(0,40000))

id <- which(grid$E == 0.40)
E40 <- AddHisto(xlist = LOSS100[id])
plot(E40, main ='', xlab ='E40 summed histograms', ylim = c(0,40000))

id <- which(grid$E == 0.75)
E75 <- AddHisto(xlist = LOSS100[id])
plot(E75, main ='', xlab ='E75 summed histograms', ylim = c(0,40000))


ratio <- grid$N/grid$Q
unique(ratio)

##### ratio #####
par(mfrow=c(2,3))
id <- which(ratio == 2)
r2 <- AddHisto(xlist = LOSS100[id])
plot(r2, main ='', xlab ='ratio 2 summed histograms', ylim = c(0,10000))

id <- which(ratio == 3)
r3 <- AddHisto(xlist = LOSS100[id])
plot(r3, main ='', xlab ='ratio 3 summed histograms', ylim = c(0,10000))

id <- which(ratio == 4)
r4 <- AddHisto(xlist = LOSS100[id])
plot(r4, main ='', xlab ='ratio 4 summed histograms', ylim = c(0,10000))


id <- which(ratio == 5)
r5 <- AddHisto(xlist = LOSS100[id])
plot(r5, main ='', xlab ='ratio 5 summed histograms', ylim = c(0,10000))

id <- which(ratio == 6)
r6 <- AddHisto(xlist = LOSS100[id])
plot(r6, main ='', xlab ='ratio 6 summed histograms', ylim = c(0,10000))


#### different scales
par(mfrow = c(1,1))
id <- which(ratio == 10)
r10 <- AddHisto(xlist = LOSS100[id])
plot(r10, main ='', xlab ='ratio 10 summed histograms', ylim = c(0,30000))

par(mfrow= c(1,2))
id <- which(ratio == 15)
r15 <- AddHisto(xlist = LOSS100[id])
plot(r15, main ='', xlab ='ratio 15 summed histograms', ylim = c(0,20000))

id <- which(ratio == 25)
r25 <- AddHisto(xlist = LOSS100[id])
plot(r25, main ='', xlab ='ratio 25 summed histograms', ylim = c(0,20000))

