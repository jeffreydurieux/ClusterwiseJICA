# Tue Dec  1 13:37:05 2020
# Author: Jeffrey Durieux, MSc


# What: summary statistics of clusterwise JICA simulation 1

load('/Volumes/LaCie/MyData/CJICA/Sim1/Sim1_results_df.Rdata')
res <- data

########### SUMMARY STATISTICS ##########

######## ARI ########
summary_by(data = res, formula = ARI ~ .)
summary_by(data = res, formula = ARI ~ Q)
summary_by(data = res, formula = ARI ~ R)
summary_by(data = res, formula = ARI ~ rho)
summary_by(data = res, formula = ARI ~ E)
summary_by(data = res, formula = ARI ~ N)

summary_by(data = res, formula = ARI ~ N:Q)
summary_by(data = res, formula = ARI ~ rho:Q)
summary_by(data = res, formula = ARI ~ E:N)

####### ARI true P #######
summary_by(data = res, formula = ARI_trueP ~ .)
summary_by(data = res, formula = ARI_trueP ~ Q)
summary_by(data = res, formula = ARI_trueP ~ R)
summary_by(data = res, formula = ARI_trueP ~ rho)
summary_by(data = res, formula = ARI_trueP ~ E)
summary_by(data = res, formula = ARI_trueP ~ N)

###################PLOTS#################

######## check aris of random100 vs trueP and perb_P #####
par(mfrow=c(1,1))
range(res$ARI_trueP)
plot(res$ARI, res$ARI_trueP, ylim = c(0,1), asp=T)
abline(0,1)

range(res$ARI_perbP)
plot(res$ARI, res$ARI_perbP, ylim = c(0,1))
abline(0,1)

######## check loss of random100 vs trueP and perb_P #####

plot(res$LOSSrandom,res$LOSStrueP)


par(mfrow= c(1,2))
LossTruep_min_Random <-  res$LOSSrandom - res$LOSStrueP 
range(LossTruep_min_Random)
plot(LossTruep_min_Random, ylim= range(LossTruep_min_Random))

LossTrueP_min_perb <- res$LOSSperb - res$LOSStrueP
plot(LossTrueP_min_perb, ylim = range(LossTruep_min_Random))

### interesting ids ####
ids <- which(LTR < 0)
res[ids,c(1,2,3,4,5, 7,8, 11,12)]


for(p in 1:100){
  cat('p index ', p, '\n')
  print(res$LOSStrueP[p])
  print(res$LOSSrandom[p])
  
}
p = 5
res$LOSStrueP[p]
res$LOSSrandom[p]



######### ANOVA ##########
library(ez)

dat <- data.frame(id = 1:4860, res)

ari_anova_res <- ezANOVA( data = dat , ARI , wid = id, within = NULL , within_full = NULL , within_covariates = NULL , between = .('Q','R','N','rho','E') , between_covariates = NULL , observed = NULL , diff = NULL , reverse_diff = FALSE , type = 3 , white.adjust = FALSE , detailed = FALSE , return_aov = TRUE)

print(ari_anova_res)
ari_anova_res$ANOVA$ges %>% round(digits = 3)
cbind(ari_anova_res$ANOVA,ef = ari_anova_res$ANOVA$ges %>% round(digits = 3))


##### ANOVA with ratio variable #####
ratio <- res$N/res$Q
res <- cbind(res,ratio)

dat <- data.frame(id = 1:4860, res)
ari_anova_res <- ezANOVA( data = dat , ARI , wid = id, within = NULL , within_full = NULL , within_covariates = NULL , between = .('ratio','E') , between_covariates = NULL , observed = NULL , diff = NULL , reverse_diff = FALSE , type = 3 , white.adjust = FALSE , detailed = FALSE , return_aov = TRUE)

print(ari_anova_res)
ari_anova_res$ANOVA$ges %>% round(digits = 3)
cbind(ari_anova_res$ANOVA,ef = ari_anova_res$ANOVA$ges %>% round(digits = 3))





##########
summary <- summary_by(data = res, formula = ARI ~ Q:N:E)
ratio <- summary$N/summary$Q
cbind(summary,ratio)

cbind(summary_by(data = res, formula = ARI ~ E:N:Q, FUN = c(mean)), ratio)
cbind(summary_by(data = res, formula = ARI ~ Q:N, FUN = c(mean,sd)), ratio)

dat <- cbind(dat, ratio = dat$N/dat$Q)
ari_anova_res <- ezANOVA( data = dat , ARI , wid = id, within = NULL , within_full = NULL , within_covariates = NULL , between = .('ratio') , between_covariates = NULL , observed = NULL , diff = NULL , reverse_diff = FALSE , type = 3 , white.adjust = FALSE , detailed = FALSE , return_aov = TRUE)

ari_anova_res

unique(dat$ratio)

library(plotly)
dat$ratio <- as.factor(dat$ratio) 
dat$E <- as.factor(dat$E) 
plot_ly(y=~ARI, x=~ratio ,type = 'box', data = dat)
plot_ly(y=~ARI, x=~E , color = ~ratio, type = 'box', data = dat) %>% 
  layout(boxmode = 'group')

plot_ly(y=~ARI, x=~ratio , color = ~E, type = 'box', data = dat) %>% 
  layout(boxmode = 'group')
