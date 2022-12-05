library(stargazer)
library(quantreg)
library(stats)
library(optim)
library(dplyr)

cps09mar <- read_dta("Documents/F22/Metrics/PS1b/cps09mar.dta")

cps09mar$exp <- cps09mar$age - 6 - cps09mar$education
data <- cps09mar[(8 < cps09mar$exp) & (cps09mar$exp < 15) & (10 < cps09mar$education) & (cps09mar$education < 18), ]
data$log_earn <- log(data$earnings)
data$exp2 <- data$exp^2
data$int <- 1

#ols regression
ols_model <- lm(log_earn ~ education + exp + exp2, data)
summary(ols_model)
stargazer(ols_model)

#quantile regression
qreg0.75 <- rq(log_earn ~ education + exp + exp2, data, tau = c(0.75))
qreg0.5 <- rq(log_earn ~ education + exp + exp2, data, tau = c(0.5))
summary(qreg0.5, se="boot")
summary(qreg0.75, se = "boot")

#gmm portion
g75_vec <- function(theta, x, W){
  fitted <-theta[1]*x$int + theta[2]*x$education + theta[3]*x$exp + theta[4]*x$exp2
  g <- (x$log_earn - fitted)*(0.75 - (x$log_earn <= fitted)) #(1-tau)*abs(x$log_earn - beta_x)* (x$log_earn <= beta_x) + tau*abs(x$log_earn - beta_x)* (x$log_earn <= beta_x)
  return(g)
}

g75 <- function(theta, x, W){
  g <- g75_vec(theta, x, W)
  m <- mean(t(g)) * W * mean(g)
  return(m)
}

g5_vec <- function(theta, x, W){
  fitted <-theta[1]*x$int + theta[2]*x$education + theta[3]*x$exp + theta[4]*x$exp2
  g <- (x$log_earn - fitted)*(0.5 - (x$log_earn <= fitted)) #(1-tau)*abs(x$log_earn - beta_x)* (x$log_earn <= beta_x) + tau*abs(x$log_earn - beta_x)* (x$log_earn <= beta_x)
  return(g)
}
g5 <- function(theta, x, W){
  g <- g5_vec(theta, x, W)
  m <- mean(t(g)) * W * mean(g)
  return(m)
}

gmm_opt75 <- optim(c(8.303, 0.154, 0.0619, -0.0017), g75, x = data, W = 1)
#ghat75 <- g75_vec(gmm_opt75$par, data, 1)
#(mean(ghat75 %*% t(ghat75)))^(-1)}
#gmm_opt_2_75 <- optim(gmm_opt75$par, g75, x = data, W = opt_w75)


gmm_opt5 <- optim(c(8.838, 0.138, -0.047, 0.0027), g5, x=data, W = 1)
#ghat5 <- g5_vec(gmm_opt5$par, data, 1)
#opt_w5 <- (mean(ghat5 %*% t(ghat5)))^(-1)
#gmm_opt_25 <- optim(gmm_opt5$par, g5, x=data, W = opt_w5)

gmm_opt75
gmm_opt5

beta <- gmm_opt75$par
gmmx <- cbind(1 + numeric(nrow(data)), data$education, data$exp, data$exp2)
omega <- 0.75*0.25*(t(gmmx) %*% gmmx)
eps <- data$log_earn - gmmx %*% beta


#use Powell estimator for Gamma (I haven't done cross validation to select the bandwidth, this is completely arbitrary)
hn <- 0.01
gamma <- numeric(4)
for (i in 1:nrow(gmmx)){
  gamma <- gamma + (abs(eps[i] < hn))*(gmmx[i,] %*% t(gmmx[i,]))/(nrow(gmmx)*hn*2)
}
#Powell method vcov
vcov <- solve(t(gamma) %*% solve(omega) %*% gamma)
vcov_powell <- xtable(vcov)
print(vcov_powell, floating=FALSE, tabular.environment="bmatrix", hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)


#standard error for beta_1 
b0_gmm_se <- sqrt(vcov[1,1]/nrow(data))
b1_gmm_se <- sqrt(vcov[2,2]/nrow(data))
b2_gmm_se <- sqrt(vcov[3,3]/nrow(data))
b3_gmm_se <- sqrt(vcov[4,4]/nrow(data))



#easy mode: here, we assume that errors are independent of regressors
f0 <- dnorm(0, mean=mean(eps), sd=sd(eps))
gamma_easy <- f0 *(t(gmmx) %*% gmmx)/sqrt(nrow(data))
vcov_easy <- solve(gamma_easy %*% solve(omega) %*% t(gamma_easy))
sqrt(vcov_easy[1,1]/nrow(data))
sqrt(vcov_easy[2,2]/nrow(data))
sqrt(vcov_easy[3,3]/nrow(data))
sqrt(vcov_easy[4,4]/nrow(data))
test <- xtable(vcov_easy)
print(test, floating=FALSE, tabular.environment="bmatrix", hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)

###Minimum distance

#data_md <- data %>% group_by(education, exp)
data$cell <- interaction(data$education, data$exp, sep="")
md_data <- data.frame(matrix(ncol = 3, nrow = length(levels(data$cell))))
colnames(md_data) <- c('log_earn_75', 'mean', 'sd')
for (i in 1:length(levels(data$cell))){
  level <- levels(data$cell)[i]
  data_sub <- subset(data, cell == level)
  fitted <- fitdistr(data_sub$log_earn, "normal")
  mu_hat <- fitted$estimate[1]
  sd_hat <- fitted$estimate[2]
  quant_75 <- quantile(data_sub$log_earn, probs=0.75)
  dens_75 <- pnorm(quant_75, mean= mu_hat, sd = sd_hat)
  new <- c(quant_75, (0.75*0.25)/(sqrt(count(data_sub)*(dens_75^2))), 0)
  md_data[i , ] <- new
}


data_md_mod <-data_md %>% group_by(cell) %>% mutate(mu = mean(log_earn - quantile(log_earn, probs=0.75)), sig = sd(log_earn - quantile(log_earn, probs=0.75)), quant_75 = quantile(log_earn, probs=0.75), n = n(), invvar =1/( 0.75*0.25/(n()*(dnorm(0, mean= mu, sd = sig))^2)))
data_md_ex <- data_md_mod[is.finite(data_md_mod$invvar),]
md_model <- lm(quant_75 ~ education + exp + exp2, data = data_md_ex, weights = invvar)
vcov_md <- vcov(md_model)


data_md2 <- aggregate(data_md_mod[, c("education", "exp", "exp2", "quant_75", "invvar")], list(data_md_mod$cell), mean)
md_model2 <- lm(quant_75 ~ education + exp + exp2, data = data_md2, weights = invvar)
vcov_md2 <- vcov(md_model2)
vcmd <- xtable(vcov_md2)
print(vcmd, floating=FALSE, tabular.environment="bmatrix", hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)
sqrt(vcov_md2[1,1]/30)
sqrt(vcov_md2[2,2]/30)
sqrt(vcov_md2[3,3]/30)
sqrt(vcov_md2[4,4]/30)



#ACTUAL MINIMUM DISTANCE CODE
mdx <- data.matrix(cbind(1 + numeric(30), data_md2[2:4]))
mdy <- data_md2$quant_75
W <- diag(data_md2$invvar)
b <- solve(t(mdx) %*% W %*% mdx) %*% t(mdx) %*% W %*% mdy
epsmd <- mdy - mdx %*% b
V <- as.numeric(var(epsmd)) * solve(t(mdx) %*% solve(W) %*% mdx)
vcmd <- xtable(V)
print(vcmd, floating=FALSE, tabular.environment="bmatrix", hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE)
sqrt(V[1,1])/30
sqrt(V[2,2])/30
sqrt(V[3,3])/30
sqrt(V[4,4])/30


N <- 400


#J = 1 estimates
sample_j <- data_md[sample(nrow(data_md), N, replace=FALSE), ]
qreg_bag <- rq(log_earn ~ education + exp + exp2, sample_j, tau = c( 0.75))
qc <- qreg_bag$coefficients[2]
gmm_75_bag <- optim(c(8.303, 0.154, 0.0619, -0.0017), g75, x = sample_j, W = 1)
gc <- gmm_75_bag$par[2]
data_md_mod <-sample_j %>% group_by(cell) %>% mutate(mu = mean(log_earn - quantile(log_earn, probs=0.75)), sig = sd(log_earn - quantile(log_earn, probs=0.75)), quant_75 = quantile(log_earn, probs=0.75), n = n(), invvar =1/( 0.75*0.25/(n()*(dnorm(0, mean= mu, sd = sig))^2)))
data_md_mod_ex <- data_md_mod[is.finite(data_md_mod$invvar),]
md_model_bag <- lm(quant_75 ~ education + exp + exp2 , data = data_md_mod_ex, weights = invvar)
mc <- md_model_bag$coefficients[2]

qc1 <- 0
gc1 <- 0
mc1 <- 0

B1 <- 1000
for (b in 1:B1){
  sample_boot <- sample_j[sample(nrow(sample_j), N, replace=TRUE), ]
  qreg_bag <- rq(log_earn ~ education + exp + exp2, sample_boot, tau = c( 0.75))
  qc1 <- qc1 + qreg_bag$coefficients[2]/B1
  gmm_75_bag <- optim(c(8.303, 0.154, 0.0619, -0.0017), g75, x = sample_boot, W = 1)
  gc1 <- gc1 + gmm_75_bag$par[2]/B1
  
  data_md_mod <-sample_boot %>% group_by(cell) %>% mutate(mu = mean(log_earn - quantile(log_earn, probs=0.75)), sig = sd(log_earn - quantile(log_earn, probs=0.75)), quant_75 = quantile(log_earn, probs=0.75), n = n(), invvar =1/( 0.75*0.25/(n()*(dnorm(0, mean= mu, sd = sig))^2)))
  data_md_mod_ex <- data_md_mod[is.finite(data_md_mod$invvar),]
  md_model_bag <- lm(quant_75 ~ education + exp + exp2 , data = data_md_mod_ex, weights = invvar)
  mc1<- mc1 + md_model_bag$coefficients[2]/B1
}


J = 400
B = 400
q_coefs <- data.frame(matrix(ncol = 1, nrow = J))
colnames(q_coefs) <- c("q_reg")
q_coefs$q_reg <- 0
q_coefs$qc <- 0
gmm_coefs <- data.frame(matrix(ncol = 1, nrow = J))
colnames(gmm_coefs) <- c("gmm_reg")
g_coefs$gc <- 0
gmm_coefs$gmm_reg<- 0
md_coefs <- data.frame(matrix(ncol = 1, nrow = J))
colnames(md_coefs) <- c("md_reg")
md_coefs$mc <- 0
md_coefs$md_reg<- 0

cumulative_time <- 0
for (j in 1:J){
  sample_j <- data_md[sample(nrow(data_md), N, replace=FALSE), ]
  old_time <- Sys.time()
  qreg_bag <- rq(log_earn ~ education + exp + exp2, sample_j, tau = c( 0.75))
  q_coefs$qc[j] <- qreg_bag$coefficients[2]
  gmm_75_bag <- optim(c(8.303, 0.154, 0.0619, -0.0017), g75, x = sample_j, W = 1)
  gmm_coefs$gc[j] <- gmm_75_bag$par[2]
  data_md_mod <-sample_j %>% group_by(cell) %>% mutate(mu = mean(log_earn - quantile(log_earn, probs=0.75)), sig = sd(log_earn - quantile(log_earn, probs=0.75)), quant_75 = quantile(log_earn, probs=0.75), n = n(), invvar =1/( 0.75*0.25/(n()*(dnorm(0, mean= mu, sd = sig))^2)))
  data_md_mod_ex <- data_md_mod[is.finite(data_md_mod$invvar),]
  md_model_bag <- lm(quant_75 ~ education + exp + exp2 , data = data_md_mod_ex, weights = invvar)
  md_coefs$mc[j] <- md_model_bag$coefficients[2]
  
  for (b in 1:B){
  sample_boot <- sample_j[sample(nrow(sample_j), N, replace=TRUE), ]
  qreg_bag <- rq(log_earn ~ education + exp + exp2, sample_boot, tau = c( 0.75))
  q_coefs$q_reg[j] <- q_coefs$q_reg[j] + qreg_bag$coefficients[2]/B
  gmm_75_bag <- optim(c(8.303, 0.154, 0.0619, -0.0017), g75, x = sample_boot, W = 1)
  #ghat75_bag <- g75_vec(gmm_75_bag$par, sample_boot, 1)
  #opt_w75_bag <- (mean(ghat75_bag %*% t(ghat75_bag)))^(-1)
  #gmm_2_75_bag <- optim(gmm_75_bag$par, g75, x = sample_boot, W = opt_w75_bag)
  gmm_coefs$gmm_reg[j] <- gmm_coefs$gmm_reg[j] + gmm_75_bag$par[2]/B
  
  data_md_mod <-sample_boot %>% group_by(cell) %>% mutate(mu = mean(log_earn - quantile(log_earn, probs=0.75)), sig = sd(log_earn - quantile(log_earn, probs=0.75)), quant_75 = quantile(log_earn, probs=0.75), n = n(), invvar =1/( 0.75*0.25/(n()*(dnorm(0, mean= mu, sd = sig))^2)))
  data_md_mod_ex <- data_md_mod[is.finite(data_md_mod$invvar),]
  md_model_bag <- lm(quant_75 ~ education + exp + exp2 , data = data_md_mod_ex, weights = invvar)
  md_coefs$md_reg[j] <- md_coefs$md_reg[j] + md_model_bag$coefficients[2]/B
  }
  new_time <- Sys.time() - old_time
  cumulative_time <-  cumulative_time + new_time
  projected <- cumulative_time/(j/J) 
  print(j/J)
  print(projected)
  print(projected - cumulative_time)
}
#started at 10:15

  
h_q <- hist(q_coefs$qc, breaks=10, xlab="Parameter estimates", main="Distribution of bagged estimates, quantile regression")
h_gmm <- hist(gmm_coefs$gc, breaks=10, xlab="Parameter estimates", main="Distribution of bagged estimates, GMM")
h_md <- hist(md_coefs$mc, breaks=20, xlab="Parameter estimates", main="Distribution of bagged estimates, minimum distance")

h_q_b <- hist(q_coefs$q_reg, breaks=10, xlab="Bootstrap estimates", main="Distribution of bootstrap estimates, quantile regression")
h_gmm_b <- hist(gmm_coefs$gmm_reg, breaks=10, xlab="Bootstrap estimates", main="Distribution of bootstrap estimates, GMM")
h_md_b <- hist(md_coefs$md_reg, breaks=10, xlab="Bootstrap estimates", main="Distribution of bootstrap estimates, minimum distance")

plot(h_q, col=rgb(1,0,0,1/4), main='Distributions of estimators', xlim=c(0.00, 0.30), ylim=c(0,130), xlab='x')
plot(h_gmm, col=rgb(0,0,1,1/4), add=TRUE)
plot(h_md, col=rgb(0,1,0,1/4),add=TRUE)
legend(0.2, 100, legend=c("qreg", "GMM", "MD"), col=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4),rgb(0,1,0,1/4)), lty=1, lwd=10)

plot(h_q_b, col=rgb(1,0,0,1/4),main='Distributions of bagged estimators', xlim=c(0.00, 0.30), ylim=c(0,110),xlab='x')
plot(h_gmm_b, col=rgb(0,0,1,1/4), add=TRUE)
plot(h_md_b, col=rgb(0,1,0,1/4),add=TRUE)
legend(0.2, 100, legend=c("qreg", "GMM", "MD"), col=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4),rgb(0,1,0,1/4)), lty=1, lwd=10)

plot(h_q, col=rgb(1,0,0,1/4),xlim=c(0.05, 0.25), ylim=c(0,110),main="Quantile regression, estimator vs bagged estimator")
plot(h_q_b, col=rgb(0,0,1, 1/4), add=TRUE)
legend(0.18, 80, legend=c("Estimator", "Bagged Estimator"), col=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), lty=1, lwd=10)

plot(h_gmm, col=rgb(1,0,0,1/4),xlim=c(0.05, 0.25), ylim=c(0,110),main="GMM, estimator vs bagged estimator")
plot(h_gmm_b, col=rgb(0,0,1, 1/4), add=TRUE)
legend(0.18, 80, legend=c("Estimator", "Bagged Estimator"), col=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), lty=1, lwd=10)

plot(h_md, col=rgb(1,0,0,1/4),xlim=c(0.05, 0.30), ylim=c(0,110),main="Minimum distance, estimator vs bagged estimator")
plot(h_md_b, col=rgb(0,0,1, 1/4), add=TRUE)
legend(0.18, 80, legend=c("Estimator", "Bagged Estimator"), col=c(rgb(1,0,0,1/4),rgb(0,0,1,1/4)), lty=1, lwd=10)



b1_q <- mean(q_coefs$qc)
b1_gmm <- mean(gmm_coefs$gc)
b1_md <- mean(md_coefs$mc)
b1_q_bag <- mean(q_coefs$q_reg)
b1_gmm_bag <- mean(gmm_coefs$gmm_reg)
b1_md_bag <- mean(md_coefs$md_reg)

true_b1 <- qreg0.75$coefficients[2]

q_coefs$bias <- q_coefs$qc - true_b1
q_coefs$bias_b <- q_coefs$q_reg - true_b1
mean(q_coefs$bias)
mean(q_coefs$bias_b)
mean(q_coefs)
sd(q_coefs$q_reg)
sd(q_coefs$qc)

gmm_coefs$bias <- gmm_coefs$gc - true_b1
gmm_coefs$bias_b <- gmm_coefs$gmm_reg - true_b1
mean(gmm_coefs$bias)
mean(gmm_coefs$bias_b)
sd(gmm_coefs$gmm_reg)
sd(gmm_coefs$gc)

md_coefs$bias <- md_coefs$mc - true_b1
md_coefs$bias_b <- md_coefs$md_reg - true_b1
mean(md_coefs$bias)
mean(md_coefs$bias_b)
sd(md_coefs$md_reg)
sd(md_coefs$mc)


