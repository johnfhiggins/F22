library(dplyr)

cereal_ps3 %>% group_by(city,quarter) %>% summarise(across(share, sum))

cereal_ps3 <- cereal_ps3 %>% 
  group_by(city,quarter) %>% 
  mutate(in_share= sum(share))

cereal_ps3$y <- log(cereal_ps3$share) - log(1-df1$in_share)

#ols_model <- lm(y ~ price + sugar + mushy, data=cereal_ps3 )
ols_model <- lm(y ~ price, data=cereal_ps3 )
summary(ols_model)

library(plm)
#ols_model_brand <- plm(y ~ price + sugar + mushy + brand , data=cereal_ps3, index=c("brand"),model="within" )
ols_model_brand <- plm(y ~ price  , data=cereal_ps3, index=c("brand"),model="within" )
summary(ols_model_brand)

library(ivreg)
inst_range <- paste("z", 1:20, sep="")
(inst_form <- as.formula(paste("y ~ price | ", paste(inst_range, collapse= "+"))))
IV_model <- ivreg(inst_form, data = cereal_ps3)
summary(IV_model)

IV_model_brand <- plm(inst_form, data=cereal_ps3, index=c("brand"), model="within")
summary(IV_model_brand)


cereal_ps3$mkup_OLS = cereal_ps3$share/(7.5856*cereal_ps3$price)
cereal_ps3$mkup_OLS_FE = cereal_ps3$share/(29.03702*cereal_ps3$price)
cereal_ps3$mkup_IV = cereal_ps3$share/(8.6859*cereal_ps3$price)
cereal_ps3$mkup_IV_FE = cereal_ps3$share/(30.1856*cereal_ps3$price)


cereal_ps3$marg_OLS = cereal_ps3$share/(7.5856)
cereal_ps3$marg_OLS_FE = cereal_ps3$share/(29.03702)
cereal_ps3$marg_IV = cereal_ps3$share/(8.6859)
cereal_ps3$marg_IV_FE = cereal_ps3$share/(30.1856)


omega = 