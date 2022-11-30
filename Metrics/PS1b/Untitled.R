library(stargazer)
library(quantreg)
library(gmm)
cps09mar <- read_dta("Documents/F22/Metrics/PS1b/cps09mar.dta")

cps09mar$exp <- cps09mar$age - 6 - cps09mar$education
data <- cps09mar[(8 < cps09mar$exp) & (cps09mar$exp < 15) & (10 < cps09mar$education) & (cps09mar$education < 18), ]
data$log_earn <- log(data$earnings)

#ols regression
ols_model <- lm(log_earn ~ education + exp + exp^2, data)
summary(ols_model)
stargazer(ols_model)

#quantile regression
qreg0.5 <- rq(log_earn ~ education + exp + exp^2, data, tau = c(0.5, 0.75))
summary(qreg0.5)

#gmm portion
