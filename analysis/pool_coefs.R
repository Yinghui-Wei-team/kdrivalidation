################################################
###              KDRI validation             ###
###   Pool coefficients across imputations   ###
################################################

#Packages
library(tidyverse)


# DCGF --------------------------------------------------------------------
#Read in data
coefs_dcgf<- read.csv("../output/recalibration/by_imputation/coefs_dcgf.csv")

#Number of imputations
n.imp<- 15

#Pool coefficients
pool_coef_dcgf<- coefs_dcgf %>%
  group_by(var) %>%
  mutate(pool_coef = mean(coef),
         var_within = sum(se.coef.^2)/n.imp,
         var_between = sum((coef - pool_coef)^2) / (n.imp-1),
         var_total = var_within + var_between + var_between/n.imp,
         pool_se = sqrt(var_total),
         lower = pool_coef - qt(0.975, df=n.imp-1)*pool_se,
         upper = pool_coef + qt(0.975, df=n.imp-1)*pool_se) %>%
  dplyr::select(-contains("var_")) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(ci = paste0("[",lower,", ",upper,"]"),
         out = paste0(pool_coef, "<br>", ci)) %>%
distinct(var, pool_coef, pool_se, ci, out)

write.csv(pool_coef_dcgf, "../output/recalibration/pool_coefs_dcgf.csv", row.names = F)



# Death as competing event ------------------------------------------------
#Read in data
coefs_csc<- read.csv("../output/recalibration/by_imputation/coefs_csc.csv")

#Number of imputations
n.imp<- 15

#Pool coefficients
pool_coef_csc<- coefs_csc %>%
  group_by(var) %>%
  mutate(pool_coef = mean(coef),
         var_within = sum(se.coef.^2)/n.imp,
         var_between = sum((coef - pool_coef)^2) / (n.imp-1),
         var_total = var_within + var_between + var_between/n.imp,
         pool_se = sqrt(var_total),
         lower = pool_coef - qt(0.975, df=n.imp-1)*pool_se,
         upper = pool_coef + qt(0.975, df=n.imp-1)*pool_se) %>%
  dplyr::select(-contains("var_")) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(ci = paste0("[",lower,", ",upper,"]"),
         out = paste0(pool_coef, "<br>", ci)) %>%
  distinct(var, pool_coef, pool_se, ci, out)

write.csv(pool_coef_csc, "../output/recalibration/pool_coefs_csc.csv", row.names = F)


