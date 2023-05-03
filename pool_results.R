################################################
###              KDRI validation             ###
###      Pool results across imputations     ###
################################################

#Packages
library(tidyverse)

#Read in output from models with results for each imputation
models<- c("cox_dcgf", "csc")
results<- c("auc_brier", "cal")
for (i in models) {
  for (j in results) {
    filepath<- file.path("../output/by_imputation", paste0(paste(i,j, sep = "_"), ".csv"))
    assign(paste(i,j,sep = "_"),
           read.csv(filepath) %>%
             mutate(model = i))
  }
}


#Collect all data frames containing calibration slope data into list
val_cal<- list()
for (i in 1:length(models)) {
  val_cal[[i]]<- get(paste(models[i],"cal", sep="_"))
  print(paste("pool",models[i],"cal", sep="_"))
  
}

val_cal<- lapply(val_cal,
                 function(x) {
                   x %>%
                     rename(est = calslope, se = cal_se) %>%
                     mutate(measure = "calslope")
                 }
)



#Collect all data frames containing AUC and Brier score data into list
val_auc_brier<- list()
for (i in 1:length(models)) {
  val_auc_brier[[i]]<- get(paste(models[i],"auc_brier", sep="_"))
  print(paste("pool",models[i],"auc_brier", sep="_"))
  
}


#Get it into correct format
val_auc_brier<- lapply(val_auc_brier,
                       function(x) {
                         x %>%
                       dplyr::select(c(imputation, model, times, AUC, Brier, se, se.1)) %>%
                       rename(AUC_est = AUC, Brier_est = Brier,
                              AUC_se = se, Brier_se = se.1) %>%
                       #Get dataframe into correct format with column for measure (AUC or Brier), estimate and se
                       pivot_longer(cols = -c(imputation, model, times),
                                    names_sep = "_",
                                    names_to = c("measure", "est_se"),
                                    values_to = "value") %>%
                       #Currently has est and se on different rows so give them thier own column
                       pivot_wider(names_from = est_se,
                                   values_from = value)
                   }
)

val_output<- reduce(c(val_cal, val_auc_brier), full_join) %>%
  group_by(model, times, measure) %>%
  summarise(pool_est = mean(est),
            var_within = sum(se^2)/n.imp,
            var_between = sum((est - pool_est)^2) / (n.imp-1),
            var_total = var_within + var_between + var_between/n.imp,
            pool_se = sqrt(var_total),
            lower = pool_est - (qnorm(0.975)*pool_se),
            upper = pool_est + (qnorm(0.975)*pool_se)) %>%
  select(-contains("var")) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(pool_est = case_when(measure=="calslope" ~ 1+pool_est,
                              TRUE ~ as.numeric(pool_est)),
         lower = case_when(measure=="calslope" ~ 1+lower,
                              TRUE ~ as.numeric(lower)),
         upper = case_when(measure=="calslope" ~ 1+upper,
                              TRUE ~ as.numeric(upper)),
         ci = paste0("[",lower,", ",upper,"]"))

write.csv(val_output, "../output/val_performance.csv", row.names = F)



################################################
###           KDRI recalibration             ###
###      Pool results across imputations     ###
################################################

#Read in output from models with results for each imputation
models<- c("cox_dcgf", "csc")
results<- c("auc", "brier", "cal", "calboot")
for (i in models) {
  for (j in results) {
    filepath<- file.path("../output/recalibration/by_imputation", paste0("rc_", paste(i,j, sep = "_"), ".csv"))
    assign(paste0("rc_", paste(i,j, sep = "_")),
           read.csv(filepath) %>%
             mutate(model = i))
  }
}

tidy_cox_cal<- rc_cox_dcgf_calboot %>%
  group_by(imputation, times) %>%
  summarise(se = (quantile(calslope,0.975) - quantile(calslope,0.025)) / (2*qnorm(0.975))) %>%
  inner_join(rc_cox_dcgf_cal %>%
              rownames_to_column(var="imputation") %>%
              pivot_longer(names_to = "times",
                           values_to = "est",
                           cols = -c(model, imputation)) %>%
              mutate(imputation = as.numeric(imputation),
                     measure = "calslope",
                     times = case_when(times=="X1" ~ 1,
                                       times=="X2" ~ 5)))
tidy_csc_cal<- rc_csc_calboot %>%
  group_by(imputation, times) %>%
  summarise(se = (quantile(calslope,0.975) - quantile(calslope,0.025)) / (2*qnorm(0.975))) %>%
  inner_join(rc_csc_cal %>%
               rownames_to_column(var="imputation") %>%
               pivot_longer(names_to = "times",
                            values_to = "est",
                            cols = -c(model, imputation)) %>%
               mutate(imputation = as.numeric(imputation),
                      measure = "calslope",
                      times = case_when(times=="X1" ~ 1,
                                        times=="X2" ~ 5)))


tidy_cox_auc<- rc_cox_dcgf_auc %>%
  rename(est = AUC) %>%
  mutate(measure = "AUC")
tidy_csc_auc<- rc_csc_auc %>%
  rename(est = AUC) %>%
  mutate(measure = "AUC")
  
tidy_cox_brier<- rc_cox_dcgf_brier %>%
  rename(est = Brier) %>%
  mutate(measure = "Brier")
tidy_csc_brier<- rc_csc_brier %>%
  rename(est = Brier) %>%
  mutate(measure = "Brier")

n.imp<- 10
recal_output_list<- list(tidy_cox_auc, tidy_cox_brier, tidy_cox_cal,
                         tidy_csc_auc, tidy_csc_brier, tidy_csc_cal)
recal_output<- reduce(recal_output_list, full_join) %>%
  group_by(model, times, measure) %>%
  summarise(pool_est = mean(est),
            var_within = sum(se^2)/n.imp,
            var_between = sum((est - pool_est)^2) / (n.imp-1),
            var_total = var_within + var_between + var_between/n.imp,
            pool_se = sqrt(var_total),
            lower = pool_est - (qnorm(0.975)*pool_se),
            upper = pool_est + (qnorm(0.975)*pool_se)) %>%
  select(-contains("var")) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(pool_est = case_when(measure=="calslope" ~ 1+pool_est,
                              TRUE ~ as.numeric(pool_est)),
         lower = case_when(measure=="calslope" ~ 1+lower,
                           TRUE ~ as.numeric(lower)),
         upper = case_when(measure=="calslope" ~ 1+upper,
                           TRUE ~ as.numeric(upper)),
         ci = paste0("[",lower,", ",upper,"]"))

write.csv(recal_output, "../output/recalibration/rc_performance.csv", row.names = F)





