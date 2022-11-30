library(stargazer)
library(quantreg)
library(stats)
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
summary(qreg0.5)
summary(qreg0.75)

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
ghat75 <- g75_vec(gmm_opt75$par, data, 1)
(mean(ghat75 %*% t(ghat75)))^(-1)}
gmm_opt_2_75 <- optim(gmm_opt75$par, g75, x = data, W = opt_w75)

benchmark("gmm1" = {gmm_opt75 <- optim(c(8.303, 0.154, 0.0619, -0.0017), g75, x = data, W = 1)},
"filler" = {ghat75 <- g75_vec(gmm_opt75$par, data, 1)
opt_w75 <- (mean(ghat75 %*% t(ghat75)))^(-1)},
"gmm2" = {gmm_opt_2_75 <- optim(gmm_opt75$par, g75, x = data, W = opt_w75)}, replications=10, columns = c("test", "replications", "elapsed",
                                                                                                            "relative", "user.self", "sys.self"))

gmm_opt5 <- optim(c(8.838, 0.138, -0.047, 0.0027), g5, x=data, W = 1)
ghat5 <- g5_vec(gmm_opt5$par, data, 1)
opt_w5 <- (mean(ghat5 %*% t(ghat5)))^(-1)
gmm_opt_25 <- optim(gmm_opt5$par, g5, x=data, W = opt_w5)


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

test <-data_md %>% group_by(cell) %>% mutate(mu = mean(log_earn), sig = sd(log_earn), quant_75 = quantile(log_earn, probs=0.75), n = n(), var = 0.75*0.25/(sqrt(n())*(dnorm(quant_75, mean= mean(log_earn), sd = sd(log_earn)))^2))
md_model_test <- lm(quant_75 ~ education + exp + exp2, data = test, weights = var)
summary(md_model_test)
N <- 400
J <- 400
B <- 400

q_coefs <- data.frame(matrix(ncol = 1, nrow = J))
colnames(q_coefs) <- c("q_reg")
q_coefs$q_reg <- 0
gmm_coefs <- data.frame(matrix(ncol = 1, nrow = J))
colnames(gmm_coefs) <- c("gmm_reg")
gmm_coefs$gmm_reg<- 0
md_coefs <- data.frame(matrix(ncol = 1, nrow = J))
colnames(md_coefs) <- c("md_reg")
md_coefs$md_reg<- 0

cumulative_time <- 0
for (j in 1:J){
  sample_j <- data_md[sample(nrow(data_md), N, replace=FALSE), ]
  old_time <- Sys.time()
  for (b in 1:B){
  sample_boot <- sample_j[sample(nrow(sample_j), N, replace=TRUE), ]
  qreg_bag <- rq(log_earn ~ education + exp + exp2, sample_boot, tau = c( 0.75))
  q_coefs[j,] <- q_coefs[j,] + qreg_bag$coefficients[2]/B
  gmm_75_bag <- optim(c(8.303, 0.154, 0.0619, -0.0017), g75, x = sample_boot, W = 1)
  #ghat75_bag <- g75_vec(gmm_75_bag$par, sample_boot, 1)
  #opt_w75_bag <- (mean(ghat75_bag %*% t(ghat75_bag)))^(-1)
  #gmm_2_75_bag <- optim(gmm_75_bag$par, g75, x = sample_boot, W = opt_w75_bag)
  gmm_coefs[j,] <- gmm_coefs[j,] + gmm_75_bag$par[2]/B
  
  data_md_mod <-sample_boot %>% group_by(cell) %>% mutate(mu = mean(log_earn - quantile(log_earn, probs=0.75)), sig = sd(log_earn - quantile(log_earn, probs=0.75)), quant_75 = quantile(log_earn, probs=0.75), n = n(), invvar =1/( 0.75*0.25/(n()*(dnorm(0, mean= mu, sd = sig))^2)))
  data_md_mod_ex <- data_md_mod[is.finite(data_md_mod$invvar),]
  md_model_bag <- lm(quant_75 ~ education + exp + exp2 , data = data_md_mod_ex, weights = invvar)
  md_coefs[j,] <- md_coefs[j,] + md_model_bag$coefficients[2]/B
  }
  new_time <- Sys.time() - old_time
  cumulative_time <-  cumulative_time + new_time
  projected <- cumulative_time/(j/J) 
  print(j/J)
  print(projected)
  print(projected - cumulative_time)
}
#started at 10:15
  
h_q <- hist(q_coefs$q_reg, breaks=10)
h_gmm <- hist(gmm_coefs$gmm_reg, breaks=10)
h_md <- hist(md_coefs$md_reg, breaks=20)

plot(h_q, col=rgb(1,0,0,1/4),main='Distribution of bootstrap estimates', xlim=c(0.08, 0.28),xlab='x')
plot(h_gmm, col=rgb(0,0,1,1/4), add=TRUE)
plot(h_md, col=rgb(0,1,0,1/4),add=TRUE)


true_b1 <- qreg0.75$coefficients[2]

q_coefs$bias <- q_coefs$q_reg - true_b1
mean(q_coefs$bias)
sd(q_coefs$q_reg)

gmm_coefs$bias <- gmm_coefs$gmm_reg - true_b1
mean(gmm_coefs$bias)
sd(gmm_coefs$gmm_reg)

md_coefs$bias <- md_coefs$md_reg - true_b1
mean(md_coefs$bias)
sd(md_coefs$md_reg)

