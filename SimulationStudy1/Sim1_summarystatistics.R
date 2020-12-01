# Tue Dec  1 13:37:05 2020
# Author: Jeffrey Durieux, MSc


# What: summary statistics of clusterwise JICA simulation 1

res <- cbind(grid, ARI = ARI)

summary_by(data = res, formula = ARI ~ .)
summary_by(data = res, formula = ARI ~ Q)
summary_by(data = res, formula = ARI ~ R)
summary_by(data = res, formula = ARI ~ rho)
summary_by(data = res, formula = ARI ~ E)

library(ez)

dat <- data.frame(id = 1:4639, res)

ari_anova_res <- ezANOVA( data = dat , ARI , wid = id, within = NULL , within_full = NULL , within_covariates = NULL , between = .('Q','R','N','rho','E') , between_covariates = NULL , observed = NULL , diff = NULL , reverse_diff = FALSE , type = 3 , white.adjust = FALSE , detailed = FALSE , return_aov = TRUE)

print(ari_anova_res)
ari_anova_res$ANOVA$ges %>% round(digits = 3)
cbind(ari_anova_res$ANOVA,ef = ari_anova_res$ANOVA$ges %>% round(digits = 3))


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
