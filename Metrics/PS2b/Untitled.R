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

cate <- data_treat %>% group_by(tr, exp) %>% mutate(n = n()) %>% group_by(tr, exp) %>% summarize(n = n(), avg=mean(log_earn))
treated <- cate[cate$tr == 1,]
untreated <- cate[cate$tr == 0,]
cond_treatment_effect <- treated$mean_log_earn - untreated$mean_log_earn 
plot(9:14, cond_treatment_effect, pch=19, xlab="Years of experience", ylab="Conditional average treatment effect")

cate_lm <- lm(log_earn ~ tr*(exp + exp2), data_treat)
summary(cate_lm)
