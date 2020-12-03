# Thu Dec  3 10:04:34 2020
# Author: Jeffrey Durieux, MSc

# summary statistics clusterwise joint ICA simulation 2

load('/Volumes/LaCie/MyData/CJICA/Sim2/Sim2_results_df.Rdata')

library(doBy)

summary_by(data = dat, f_SpecKernARI ~ ., FUN = c(mean,sd), na.rm = T)
summary_by(data = dat, f_PamKernARI ~ ., FUN = c(mean,sd), na.rm = T)
summary_by(data = dat, f_HclKernARI ~ ., FUN = c(mean,sd), na.rm = T)
summary_by(data = dat, f_PamCosARI ~ ., FUN = c(mean,sd), na.rm = F)
summary_by(data = dat, f_HclCosARI ~ ., FUN = c(mean,sd), na.rm = F)

summary_by(data = dat, icaQ_SpecKernARI ~ ., FUN = c(mean,sd), na.rm = F)
summary_by(data = dat, icaQ_PamKernARI ~ ., FUN = c(mean,sd), na.rm = F)
summary_by(data = dat, icaQ_HclKernARI ~ ., FUN = c(mean,sd), na.rm = F)
summary_by(data = dat, icaQ_PamCosARI ~ ., FUN = c(mean,sd), na.rm = F)
summary_by(data = dat, icaQ_HclCosARI ~ ., FUN = c(mean,sd), na.rm = F)

summary_by(data = dat, icaQR_SpecKernARI ~ ., FUN = c(mean,sd), na.rm = F)
summary_by(data = dat, icaQR_PamKernARI ~ ., FUN = c(mean,sd), na.rm = F)
summary_by(data = dat, icaQR_HclKernARI ~ ., FUN = c(mean,sd), na.rm = F)
summary_by(data = dat, icaQR_PamCosARI ~ ., FUN = c(mean,sd), na.rm = F)
summary_by(data = dat, icaQR_HclCosARI ~ ., FUN = c(mean,sd), na.rm = F)


load('/Volumes/LaCie/MyData/CJICA/Sim1/Sim1_results_df.Rdata')

par(mfrow = c(1,3))

#### FULL data ####

# spectral clustering rbf
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$f_SpecKernLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff,ylim = c(0,2500))

dd <- data$LOSSrandom - dat$f_SpecKernLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$f_SpecKernLOSS, asp = F, col = col,main = table(ss));abline(0,1)

#pam rbf
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$f_PamKernLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff, ylim = c(0,2500))

dd <- data$LOSSrandom - dat$f_PamKernLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$f_PamKernLOSS, asp = F, col = col,main = table(ss));abline(0,1)

# hcl rbf
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$f_HclKernLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff, ylim = c(0,2500))

dd <- data$LOSSrandom - dat$f_HclKernLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$f_HclKernLOSS, asp = F, col = col,main = table(ss));abline(0,1)

# pam cos
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$f_PamCosLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff, ylim = c(0,2500))

dd <- data$LOSSrandom - dat$f_PamCosLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$f_PamCosLOSS, asp = F, col = col,main = table(ss));abline(0,1)

# hcl cos
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$f_HclCosLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff, ylim = c(0,2500))

dd <- data$LOSSrandom - dat$f_HclCosLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$f_HclCosLOSS, asp = F, col = col,main = table(ss));abline(0,1)

######### ICA Q ###########

# spectral clustering rbf
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$icaQ_SpecKernLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff,ylim = c(0,2500))

dd <- data$LOSSrandom - dat$icaQ_SpecKernLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$icaQ_SpecKernLOSS, asp = F, col = col,main = table(ss));abline(0,1)

#pam rbf
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$icaQ_PamKernLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff, ylim = c(0,2500))

dd <- data$LOSSrandom - dat$icaQ_PamKernLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$icaQ_PamKernLOSS, asp = F, col = col,main = table(ss));abline(0,1)

# hcl rbf
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$icaQ_HclKernLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff, ylim = c(0,2500))

dd <- data$LOSSrandom - dat$icaQ_HclKernLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$icaQ_HclKernLOSS, asp = F, col = col,main = table(ss));abline(0,1)

# pam cos
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$icaQ_PamCosLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff, ylim = c(0,2500))

dd <- data$LOSSrandom - dat$icaQ_PamCosLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$icaQ_PamCosLOSS, asp = F, col = col,main = table(ss));abline(0,1)

# hcl cos
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$icaQ_HclCosLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff, ylim = c(0,2500))

dd <- data$LOSSrandom - dat$icaQ_HclCosLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$icaQ_HclCosLOSS, asp = F, col = col,main = table(ss));abline(0,1)


######### ICA QR ###########

# spectral clustering rbf
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$icaQR_SpecKernLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff,ylim = c(0,2500))

dd <- data$LOSSrandom - dat$icaQR_SpecKernLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$icaQR_SpecKernLOSS, asp = F, col = col,main = table(ss));abline(0,1)

#pam rbf
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$icaQR_PamKernLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff, ylim = c(0,2500))

dd <- data$LOSSrandom - dat$icaQR_PamKernLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$icaQR_PamKernLOSS, asp = F, col = col,main = table(ss));abline(0,1)

# hcl rbf
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$icaQR_HclKernLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff, ylim = c(0,2500))

dd <- data$LOSSrandom - dat$icaQR_HclKernLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$icaQR_HclKernLOSS, asp = F, col = col,main = table(ss));abline(0,1)

# pam cos
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$icaQR_PamCosLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff, ylim = c(0,2500))

dd <- data$LOSSrandom - dat$icaQR_PamCosLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$icaQR_PamCosLOSS, asp = F, col = col,main = table(ss));abline(0,1)

# hcl cos
LOSSdiff_random <- data$LOSSrandom - data$LOSStrueP
LOSSdiff <- dat$icaQR_HclCosLOSS - data$LOSStrueP
plot(LOSSdiff_random, ylim = c(0,2500))
plot(LOSSdiff, ylim = c(0,2500))

dd <- data$LOSSrandom - dat$icaQR_HclCosLOSS
ss <- sign(dd)
table(ss)
col <- ifelse(dd > 0, 'red', 'darkgreen')
plot(data$LOSSrandom, dat$icaQR_HclCosLOSS, asp = F, col = col,main = table(ss));abline(0,1)



###### LR on diff ######
#1 = rational is equal to random or winner
#0 = random is the winner

dd <- data$LOSSrandom - dat$icaQR_PamCosLOSS
ss <- sign(dd)
table(ss)
idd <- ss==0
ss[idd] <- 1
idd <- ss==-1
ss[idd] <- 0
table(ss)
#ss <- ifelse(test = ss == 0, 1, 0)

df <- cbind(dat[,1:5],ss)

library(broom)
library(plotly)
LR <- glm(data = df, formula = ss ~ ., family = 'binomial') 
LR
tidy(LR)
exp(coef(LR))

#2nd order int
LR2 <- glm(data = df, formula = ss ~ . +.^2, family = 'binomial')

LR2
summary(LR2)
tidy(LR2)
exp(coef(LR2))

anova(LR,LR2)

#3nd order int
LR3 <- glm(data = df, formula = ss ~ . + .^3, family = 'binomial') 
LR3
summary(LR3)
tidy(LR3)
exp(coef(LR3)) %>%round(digits = 3)

anova(LR,LR2)


#### 
ratio <- df$N/df$Q
df <- cbind(df, ratio)
LRratio <- glm(data = df, formula = ss ~ ratio, family = 'binomial')
summary(LRratio)
tidy(LRratio)
exp(coef(LRratio)) %>%round(digits = 3)
