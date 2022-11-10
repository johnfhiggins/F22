library(stats)
library(MASS)
library(stargazer)

Xmat <- read_csv("Xmat.csv")
probit_wm <- glm(walmart ~ log_pop + log_ret + urban_pct + midwest + south + kmart, family = binomial(link = "probit"), data = Xmat)
probit_km <- glm(kmart ~ log_pop + log_ret + urban_pct + midwest + south + walmart, family = binomial(link = "probit"), data = Xmat)
summary(probit_wm)
summary(probit_km)
stargazer(probit_wm, probit_km)

#pred_wm = predict(probit_wm, Xmat, type="response")
#pred_km = predict(probit_km, Xmat, type="response")
#hist(as.numeric(pred_wm))
#hist(as.numeric(pred_km))
#Xmat$pred_wm <- as.numeric(pred_wm)
#Xmat$pred_km <- as.numeric(pred_km)


Xmat$large_stores <- Xmat$walmart + Xmat$kmart
Xmat$all_stores <- Xmat$large_stores + Xmat$small_stores
Xmat$large_st_ordered <- as.ordered(Xmat$large_stores)
Xmat$all_st_ordered <- as.ordered(Xmat$all_stores)
ordered_prob_WK <- polr(large_st_ordered ~ log_pop +  log_ret + urban_pct + midwest + south, data=Xmat, method= c("probit"))
ordered_prob_all <- polr(all_st_ordered ~ log_pop +  log_ret + urban_pct + midwest + south, data=Xmat, method= c("probit"))
summary(ordered_prob_WK)
summary(ordered_prob_all)

#Xmat$pih_w <- qnorm(Xmat$pred_wm) #log(Xmat$pred_wm) - log(1-Xmat$pred_wm)
#Xmat$pih_k <- qnorm(Xmat$pred_km) #log(Xmat$pred_km) - log(1-Xmat$pred_km)

#walmart_reg <- lm(pih_w ~ log_pop + log_ret + urban_pct + midwest + south+ kmart, Xmat)
#kmart_reg <- lm(pih_k ~ log_pop + log_ret + urban_pct + midwest + south+ walmart, Xmat)
#summary(walmart_reg)
#summary(kmart_reg)


#probit with benton county distance instrument
#w_reg_benton <- glm(walmart ~ log_pop + log_ret + log_benton_dist + urban_pct + midwest + south, family=binomial(link="probit"), data=Xmat)
#pred_wm_bent = predict(w_reg_benton, Xmat, type="response")
#Xmat$pred_wm_bent <- as.numeric(pred_wm_bent)

#first stage regression 
w_reg <- glm(walmart ~ log_pop  + log_ret + log_benton_dist + urban_pct +  midwest + south  , family = binomial(link = "probit"), data = Xmat)
k_reg <- glm(kmart ~ log_pop  + log_ret + urban_pct + midwest + south , family = binomial(link = "probit"), data = Xmat)
summary(w_reg)
summary(k_reg)



#find probit predictions
pred_wm = predict(w_reg, Xmat, type="response")
pred_km = predict(k_reg, Xmat, type="response")
hist(as.numeric(pred_wm))
hist(as.numeric(pred_km))
#add to dataframe
Xmat$pred_wm <- as.numeric(pred_wm)
Xmat$pred_km <- as.numeric(pred_km)
#implied profit functions based on predicted entry choice
Xmat$pih_w <- qnorm(Xmat$pred_wm)
Xmat$pih_k <- qnorm(Xmat$pred_km) 


#regress implied profits on county covariates and the predicted probability of entry for the rival. 
wm_step2_prob <- lm(pih_w ~ log_pop  + log_ret +log_benton_dist+ urban_pct + midwest + south  + pred_km, data = Xmat)
km_step2_prob <- lm(pih_k ~ log_pop  + log_ret + urban_pct + midwest + south + pred_wm, data = Xmat)
summary(wm_step2_prob)
summary(km_step2_prob)
#using the above regressions, find the WLS weight matrix and do weighted least squares
wt <- 1 / lm(abs(wm_step2_prob$residuals) ~ wm_step2_prob$fitted.values)$fitted.values^2
wls_prob_w <- lm(pih_w ~ log_pop  + log_ret + log_benton_dist + urban_pct + midwest + south  + pred_km, data = Xmat, weights = wt)
summary(wls_prob_w)
wt2 <- 1 / lm(abs(km_step2_prob$residuals) ~ km_step2_prob$fitted.values)$fitted.values^2
wls_prob_k <- lm(pih_k ~ log_pop  + log_ret + urban_pct + midwest + south  + pred_wm, data = Xmat, weights = wt2)
summary(wls_prob_k)
stargazer(wls_prob_w, wls_prob_k)


#create moment functions

g1 <- function(theta, x){
  m1 <- x$walmart - sigma_w(theta, x)
  return(m1)
}
g2 <- function(theta, x){
  m1 <- x$kmart - sigma_k(theta, x)
  return(m1)
}
sigma_w <- function(theta, x){
  pi <- theta[1] + x$log_pop*theta[2] + x$log_ret*theta[3] + x$log_benton_dist*theta[4] + x$urban_pct*theta[5] + x$midwest*theta[6] + x$south*theta[7] + x$pred_km*theta[8]
  prob <- pnorm(pi)
  return(prob)
}
sigma_k <- function(theta, x){
  pi <- theta[1] + x$log_pop*theta[2] + x$log_ret*theta[3]  + x$urban_pct*theta[4] + x$midwest*theta[5] + x$south*theta[6] + x$pred_wm*theta[7]
  prob <- pnorm(pi)
  return(prob)
}

#gmm computation
gmm_wm_2 <- gmm(g1, x=Xmat, t0 = c(-12.2, 1.823,1.544,-1.079,1.306,0.0,0.702,-0.714))
gmm_km_2 <- gmm(g2, x=Xmat, t0 = c(-21.6, 1.758,1.730,1.236,0.58,0.24,-0.76))
summary(gmm_wm_2)
summary(gmm_km_2)
stargazer(gmm_wm_2)


#fully parametric specification
prob_wm <- glm(walmart ~ log_pop + log_ret + log_benton_dist + urban_pct + midwest + south + pred_km, data=Xmat, family = binomial(link = "probit"))
prob_km <- glm(kmart ~ log_pop + log_ret  + urban_pct + midwest + south + pred_wm, data=Xmat, family = binomial(link = "probit"))
summary(prob_wm)
summary(prob_km)


##Ignore all of the stuff below here! I attempted to do nonparemetric estimation for the entry probability, but it didn't quite work out too well

wm_gmm <- gmm(x_dat$pih_w ~ x_dat$log_pop + x_dat$log_ret +x_dat$log_benton_dist+ x_dat$urban_pct + x_dat$midwest + x_dat$south + x_dat$pred_wm, x=x_dat)
km_gmm <- gmm(x_dat$pih_w ~ x_dat$log_pop + x_dat$log_ret + x_dat$urban_pct + x_dat$midwest + x_dat$south + x_dat$pred_km, x=x_dat)
summary(wm_gmm)
summary(km_gmm)

library(Sieve)
sieve.model_wm <- sieve_preprocess(X = Xmat[, c('log_pop', 'log_ret', 'log_benton_dist',  'urban_pct', 'midwest', 'south')], basisN = 500, type='cosine')
sieve.fit_wm <- sieve_solver(model = sieve.model_wm, Y = Xmat$walmart, family='binomial')
s_pred_wm <- sieve_predict(model=sieve.fit_wm, testX = Xmat[, c('log_pop', 'log_ret', 'log_benton_dist',  'urban_pct', 'midwest', 'south')])
wm_sieve <- unlist(data.frame(s_pred_wm[1])[93])

sieve.model_km <- sieve_preprocess(X = Xmat[, c('log_pop', 'log_ret',  'urban_pct', 'midwest', 'south')], basisN = 500, type='cosine')
sieve.fit_km <- sieve_solver(model = sieve.model_km, Y = Xmat$kmart, family='binomial')
s_pred_km <- sieve_predict(model=sieve.fit_km, testX = Xmat[, c('log_pop', 'log_ret',  'urban_pct', 'midwest', 'south')])
km_sieve <- unlist(data.frame(s_pred_km[1])[100])

Xmat$sieve_wm <- wm_sieve
Xmat$sieve_km <- km_sieve
Xmat$pih_w_s <- qnorm(Xmat$sieve_wm)
Xmat$pih_k_s <- qnorm(Xmat$sieve_km)
wm_sieve_reg <- lm(pih_w_s ~ log_pop  + log_ret + log_benton_dist + urban_pct + midwest + south  + sieve_km, data = Xmat)
km_sieve_reg <- lm(pih_k_s ~ log_pop  + log_ret + urban_pct + midwest + south + sieve_wm, data = Xmat)
wtws <- 1 / lm(abs(wm_sieve_reg$residuals) ~ wm_sieve_reg$fitted.values)$fitted.values^2
wls_sieve_w <- lm(pih_w_s ~ log_pop  + log_ret + log_benton_dist + urban_pct + midwest + south  + sieve_km, data = Xmat, weights = wtws)
summary(wls_sieve_w)
wtks <- 1 / lm(abs(km_sieve_reg$residuals) ~ km_sieve_reg$fitted.values)$fitted.values^2
wls_sieve_k<- lm(pih_k_s ~ log_pop  + log_ret  + urban_pct + midwest + south  + sieve_wm, data = Xmat, weights = wtks)
summary(wls_sieve_k)

summary(wm_sieve_reg)
summary(km_sieve_reg)
stargazer(wm_step2_reg, wm_sieve_reg, km_step2_reg, km_sieve_reg)

library(gmm)
x_dat <- data.frame(Xmat[,c('log_pop', 'log_ret','log_benton_dist', 'urban_pct', 'midwest', 'south', 'sieve_wm', 'sieve_km', 'pred_km','pred_wm', 'pih_w_s', 'pih_k_s')])

gmm_wm <- gmm(x_dat$pih_w_s ~ x_dat$log_pop + x_dat$log_ret +x_dat$log_benton_dist+ x_dat$urban_pct + x_dat$midwest + x_dat$south + x_dat$sieve_km, x=x_dat, type="iterative")
gmm_km <- gmm(x_dat$pih_k_s ~ x_dat$log_pop + x_dat$log_ret + x_dat$urban_pct + x_dat$midwest + x_dat$south + x_dat$sieve_wm, x=x_dat, type="iterative")
summary(gmm_wm)
summary(gmm_km)

#also do this for logit just for fun
Xmat$pih_w_s_log <- log(Xmat$sieve_wm) - log(1-Xmat$sieve_wm)
Xmat$pih_k_s_log <- log(Xmat$sieve_km) - log(1 - Xmat$sieve_km)

library(gmm)
x_dat <- data.frame(Xmat[,c('log_pop', 'log_ret','log_benton_dist', 'urban_pct', 'midwest', 'south', 'sieve_wm', 'sieve_km', 'pred_km','pred_wm', 'pih_w_s', 'pih_k_s', 'pih_w_s_log', 'pih_k_s_log')])

gmm_wm_log <- gmm(x_dat$pih_w_s_log ~ x_dat$log_pop + x_dat$log_ret +x_dat$log_benton_dist+ x_dat$urban_pct + x_dat$midwest + x_dat$south + x_dat$sieve_km, x=x_dat)
gmm_km_log <- gmm(x_dat$pih_k_s_log ~ x_dat$log_pop + x_dat$log_ret + x_dat$urban_pct + x_dat$midwest + x_dat$south + x_dat$sieve_wm, x=x_dat)
summary(gmm_wm_log)
summary(gmm_km_log)
stargazer(gmm_wm, gmm_wm_log, gmm_km, gmm_km_log)

