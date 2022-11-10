library(stats)
library(MASS)
library(stargazer)

Xmat <- read_csv("Xmat.csv")
probit_wm <- glm(walmart ~ log_pop  + urban_pct + midwest + south + kmart, family = binomial(link = "probit"), data = Xmat)
probit_km <- glm(kmart ~ log_pop  + urban_pct + midwest + south + walmart, family = binomial(link = "probit"), data = Xmat)
summary(probit_wm)
summary(probit_km)
stargazer(probit_wm, probit_km)

pred_wm = predict(probit_wm, Xmat, type="response")
pred_km = predict(probit_km, Xmat, type="response")
hist(as.numeric(pred_wm))
hist(as.numeric(pred_km))
Xmat$pred_wm <- as.numeric(pred_wm)
Xmat$pred_km <- as.numeric(pred_km)


Xmat$large_stores <- Xmat$walmart + Xmat$kmart
Xmat$all_stores <- Xmat$large_stores + Xmat$small_stores
Xmat$large_st_ordered <- as.ordered(Xmat$large_stores)
Xmat$all_st_ordered <- as.ordered(Xmat$all_stores)
ordered_prob_WK <- polr(large_st_ordered ~ log_pop  + urban_pct + midwest + south, data=Xmat, method= c("probit"))
ordered_prob_all <- polr(all_st_ordered ~ log_pop  + urban_pct + midwest + south, data=Xmat, method= c("probit"))
summary(ordered_prob_WK)
summary(ordered_prob_all)

Xmat$pih_w <- qnorm(Xmat$pred_wm) #log(Xmat$pred_wm) - log(1-Xmat$pred_wm)
Xmat$pih_k <- qnorm(Xmat$pred_km) #log(Xmat$pred_km) - log(1-Xmat$pred_km)

walmart_reg <- lm(pih_w ~ log_pop  + urban_pct + midwest + south+ pred_km, Xmat)
kmart_reg <- lm(pih_k ~ log_pop  + urban_pct + midwest + south+ pred_wm, Xmat)
summary(walmart_reg)
summary(kmart_reg)



#probit with benton county distance IV ???

w_reg <- glm(walmart ~ log_pop   + log_benton_dist + urban_pct + midwest + south , family = binomial(link = "probit"), data = Xmat)
k_reg <- glm(kmart ~ log_pop   + urban_pct + midwest + south, family = binomial(link = "probit"), data = Xmat)
#library(ivprobit)
#iv_prob_KM <- ivprobit(kmart ~ log_pop + log_ret + urban_pct + midwest + south | walmart | log_pop + log_ret + urban_pct + midwest + south +log_benton_dist, Xmat)
##summary(iv_prob_KM)
#pred_km_prob <- fitted(iv_prob_KM, Xmat, type="response")
pred_wm = predict(w_reg, Xmat, type="response")
pred_km = predict(k_reg, Xmat, type="response")
hist(as.numeric(pred_wm))
hist(as.numeric(pred_km))
Xmat$pred_wm <- as.numeric(pred_wm)
Xmat$pred_km <- as.numeric(pred_km)
Xmat$pih_w <- qnorm(Xmat$pred_wm) #log(Xmat$pred_wm) - log(1-Xmat$pred_wm)
Xmat$pih_k <- qnorm(Xmat$pred_km) 

wm_step2_reg <- lm(pih_w ~ log_pop  + log_ret + log_benton_dist + urban_pct + midwest + south  + pred_km, data = x_dat)
km_step2_reg <- lm(pih_k ~ log_pop  + log_ret + urban_pct + midwest + south + pred_wm, data = x_dat)
summary(wm_step2_reg)
summary(km_step2_reg)

library(Sieve)
sieve.model_wm <- sieve_preprocess(X = Xmat[, c('log_pop',  'log_benton_dist',  'urban_pct', 'midwest', 'south')], basisN = 50, type='cosine')
sieve.fit_wm <- sieve_solver(model = sieve.model_wm, Y = Xmat$walmart, family='binomial')
s_pred_wm <- sieve_predict(model=sieve.fit_wm, testX = Xmat[, c('log_pop', 'log_benton_dist',  'urban_pct', 'midwest', 'south')])
wm_sieve <- unlist(data.frame(s_pred_wm[1])[86])

sieve.model_km <- sieve_preprocess(X = Xmat[, c('log_pop',  'urban_pct', 'midwest', 'south')], basisN = 50, type='cosine')
sieve.fit_km <- sieve_solver(model = sieve.model_km, Y = Xmat$kmart, family='binomial')
s_pred_km <- sieve_predict(model=sieve.fit_km, testX = Xmat[, c('log_pop',  'urban_pct', 'midwest', 'south')])
km_sieve <- unlist(data.frame(s_pred_km[1])[90])

Xmat$sieve_wm <- wm_sieve
Xmat$sieve_km <- km_sieve
Xmat$pih_w_s <- qnorm(Xmat$sieve_wm)
Xmat$pih_k_s <- qnorm(Xmat$sieve_km)
wm_sieve_reg <- lm(pih_w_s ~ log_pop   + log_benton_dist + urban_pct + midwest + south  + sieve_km, data = Xmat)
km_sieve_reg <- lm(pih_k_s ~ log_pop   + urban_pct + midwest + south + sieve_wm, data = Xmat)
summary(wm_sieve_reg)
summary(km_sieve_reg)
stargazer(wm_step2_reg, wm_sieve_reg, km_step2_reg, km_sieve_reg)

library(gmm)
x_dat <- data.frame(Xmat[,c('log_pop', 'log_benton_dist', 'urban_pct', 'midwest', 'south', 'sieve_wm', 'sieve_km','kmart', 'walmart', 'pred_km','pred_wm','pih_w', 'pih_k', 'pih_w_s', 'pih_k_s')])
#x_dat_km <- data.frame(Xmat[,c('log_pop', 'log_ret', 'urban_pct', 'midwest', 'south', 'sieve_km', 'pred_wm', 'pih_k_s')])
moments <- function(theta, data) {
  #thet <- as.matrix(theta)
  y <- as.matrix(data[, 7])
  x <- data[, 1:6]
  m <- (y -x%*%theta) 
  return(cbind(m))
}
gmm_wm <- gmm(x_dat$pih_w_s ~ x_dat$log_pop  +x_dat$log_benton_dist+ x_dat$urban_pct + x_dat$midwest + x_dat$south + x_dat$sieve_km, x=x_dat, type="iterative", prewhite=0)
gmm_km <- gmm(x_dat$pih_k_s ~ x_dat$log_pop + x_dat$urban_pct + x_dat$midwest + x_dat$south + x_dat$sieve_wm, x=x_dat, type="iterative", prewhite=0)
summary(gmm_wm)
summary(gmm_km)
stargazer(wm_step2_reg, gmm_wm, km_step2_reg, gmm_km)


