library(stargazer)
library(quantreg)
library(stats)
library(optim)
library(dplyr) 
library(haven)

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

data_treat <- data[(data$education == 12 | data$education == 16),]
data_treat$tr <- as.numeric(data_treat$education ==16)

#full sample computation

#naive ate
treated <- data_treat[data_treat$tr ==1,]
untreated <- data_treat[data_treat$tr==0,]
ate <- mean(treated$log_earn) - mean(untreated$log_earn)

#treatment effect regression
treat_lm <- lm(log_earn ~ tr*(exp + exp2), data=data_treat)
ate_reg <- treat_lm$coefficients[2] + treat_lm$coefficients[5]*mean(data_treat$exp) + treat_lm$coefficients[6]*mean(data_treat$exp2)

#propensity matching

prop_data <- data_treat %>% group_by(tr, exp) %>% mutate(n = n()) %>% group_by(tr, exp) %>% summarize(n = n(), mean_earn=mean(log_earn))
prop_score <- prop_data %>% group_by(exp) %>% mutate(P = n/sum(n), nexp =sum(n)) %>% group_by(tr,exp) %>% summarize(P, nexp, avg=mean(mean_earn))
prop_score$P1 <- ifelse(prop_score$tr==1, prop_score$P, 1-prop_score$P)
prop_ordered <- prop_score[order(prop_score$P1),]
treated <- prop_ordered[prop_ordered$tr==1,]
untreated <- prop_ordered[prop_ordered$tr==0,]
exp_lower <- treated$exp[1:2]
exp_middle <- treated$exp[3:4]
exp_high <- treated$exp[5:6]
data_lower <- data_treat[data_treat$exp %in% exp_lower,]
data_middle <- data_treat[data_treat$exp %in% exp_middle,]
data_high <- data_treat[data_treat$exp %in% exp_high,]
ate_lower <- ate_naive(data_lower)
ate_middle <- ate_naive(data_middle)
ate_high <- ate_naive(data_high)
ate_prop <- mean(ate_lower, ate_middle, ate_high)

ate_naive <- function(data_naive){
  treated <- data_naive[data_naive$tr ==1,]
  untreated <- data_naive[data_naive$tr==0,]
  ate <- mean(treated$log_earn) - mean(untreated$log_earn)
  ate
}

ate_reg <- function(data_reg){
  treat_lm <- lm(log_earn ~ tr*(exp + exp2), data=data_reg)
  ate_reg <- treat_lm$coefficients[2] + treat_lm$coefficients[5]*mean(data_reg$exp) + treat_lm$coefficients[6]*mean(data_reg$exp2)
}

ate_prop <- function(data_prop){
  prop_data <- data_prop %>% group_by(tr, exp) %>% mutate(n = n()) %>% group_by(tr, exp) %>% summarize(n = n(), mean_earn=mean(log_earn))
  prop_score <- prop_data %>% group_by(exp) %>% mutate(P = n/sum(n), nexp =sum(n)) %>% group_by(tr,exp) %>% summarize(P, nexp, avg=mean(mean_earn))
  prop_score$P1 <- ifelse(prop_score$tr==1, prop_score$P, 1-prop_score$P)
  prop_ordered <- prop_score[order(prop_score$P1),]
  treated <- prop_ordered[prop_ordered$tr==1,]
  untreated <- prop_ordered[prop_ordered$tr==0,]
  exp_lower <- treated$exp[1:2]
  exp_middle <- treated$exp[3:4]
  exp_high <- treated$exp[5:6]
  data_lower <- data_prop[data_prop$exp %in% exp_lower,]
  data_middle <- data_prop[data_prop$exp %in% exp_middle,]
  data_high <- data_prop[data_prop$exp %in% exp_high,]
  ate_lower <- ate_naive(data_lower)
  ate_middle <- ate_naive(data_middle)
  ate_high <- ate_naive(data_high)
  ate_mean <- mean(ate_lower, ate_middle, ate_high)
  return(ate_mean)
}

#gets rid of annoying warning
options(dplyr.summarise.inform = FALSE)


N <- 400
B <- 500

#generate 500 random samples of length 400 to be used by each estimation method. We want to use the same draws across our three samples to avoid simulation error impacting our comparisons of different estimators
data_400 <-data_treat[(sample(nrow(data_treat), N, replace=FALSE)),]

ate_n_400 <- ate_naive(data_400)
ate_r_400 <- ate_reg(data_400)
ate_p_400 <- ate_prop(data_400)

sample_indices <- replicate(B, sample(nrow(data_treat), N, replace=FALSE))
ate_n <- numeric(B)
ate_r <- numeric(B)
ate_p <- numeric(B)
for (i in 1:B){
  sample <- data_treat[sample_indices[i,],]
  ate_n[i] <- ate_naive(sample)
  ate_r[i] <- ate_reg(sample)
  ate_p[i] <- ate_prop(sample)
}
hist_df <- data.frame(var = c(rep('Naive', B), rep('Regression', B), rep('Propensity', B)), value=c(ate_n, ate_r, ate_p))

ggplot(subset(hist_df, hist_df$var=='Naive'), aes(x = value, fill=var)) + geom_histogram(fill='blue', alpha=0.5, position='identity') + xlab("Estimated ATE") + labs(title = "Histogram of ATE estimates, Naive ATE", subtitle="Estimates based on 500 samples of size 400 (vertical lines are full sample OLS and ATE estimates)") + geom_vline(xintercept=0.6037) + geom_vline(xintercept=0.6055)

ggplot(subset(hist_df, hist_df$var=='Regression'), aes(x = value, fill=var)) + geom_histogram(fill='blue', alpha=0.5, position='identity') + xlab("Estimated ATE") + labs(title = "Histogram of ATE estimates, treatment effects regression", subtitle="Estimates based on 500 samples of size 400 (vertical lines are full sample OLS and ATE estimates)") + geom_vline(xintercept=0.6037) + geom_vline(xintercept=0.6055)
ggplot(subset(hist_df, hist_df$var=='Propensity'), aes(x = value, fill=var)) + geom_histogram(fill='blue', alpha=0.5, position='identity') + xlab("Estimated ATE") + labs(title = "Histogram of ATE estimates, propensity score", subtitle="Estimates based on 500 samples of size 400 (vertical lines are full sample OLS and ATE estimates)") + geom_vline(xintercept=0.6037) + geom_vline(xintercept=0.6055)

ggplot(hist_df, aes(x = value, fill=var)) + geom_histogram(color='#e9ecef', alpha=0.45, position='identity') + xlab("Estimated ATE") + labs(title = "Histogram of ATE estimates by estimation method", subtitle="Estimates based on 500 samples of size 400", fill = "Estimation method")

mean(ate_n)
sd(ate_n)
mean(ate_r)
sd(ate_r)
mean(ate_p)
sd(ate_p)

