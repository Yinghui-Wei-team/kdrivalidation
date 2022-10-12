################################################
###            KDRI recalibration            ###
###                 Cox model                ###
################################################

#Packages
library(tidyverse)
library(survival)
library(riskRegression)
library(pec)
library(prodlim)

#Data
data<- read.csv("../../NHSBT Data/Data/kdri_mi.csv", na.strings = "")

D<- data %>%
  filter(X_mi_m > 0 & X_mi_m <=10) %>%
  mutate(gsurv = case_when(gsurv==0 ~ 0.5,
                           TRUE ~ as.numeric(gsurv)),
         gsurv_yr=gsurv/365.5) %>%
  select(arecip_id, gsurv_yr, gcens, uskdri_lp, X_mi_m,
         dage, dage18_i, dage50_i, dheight, dweight, dweight80_i,
         dethnicb_i, dpast_hypertension, dpast_diabetes, cva_i,
         dcreatmgdl, dcreat15_i, dhcv, dtype) %>%
  mutate(rcdage18 = (dage-18)*dage18_i, rcdage50 = (dage-50)*dage50_i, rcdage = dage-40,
         rcdheight = (dheight-170)/10,
         rcdweight = (dweight-80)/5*dweight80_i,
         rcdcreat = dcreatmgdl-1, rcdcreat15 = (dcreatmgdl-1.5)*dcreat15_i,
         #add indicator function for calculating linear predictor
         dpast_hypertension_i = case_when(dpast_hypertension == "Yes" ~ 1,
                                          dpast_hypertension == "No" ~ 0),
         dpast_diabetes_i = case_when(dpast_diabetes == "Yes" ~ 1,
                                      dpast_diabetes == "No" ~ 0),
         dhcv_i = case_when(dhcv == "Positive" ~ 1,
                            dhcv == "Negative" ~ 0),
         dtype_i = case_when(dtype == "DCD" ~ 1,
                             dtype == "DBD" ~ 0))

#Recalibrate the KDRI using the Cox model
#Number of imputations
n.imp<- 10

#Model formula
outcome<- "Surv(gsurv_yr, gcens)"
vars<- c("rcdage18", "rcdage", "rcdage50", "rcdheight", "rcdweight",
         "dethnicb_i", "dpast_hypertension_i", "dpast_diabetes_i", "cva_i",
         "rcdcreat", "rcdcreat15", "dhcv_i", "dtype_i")

kdriform<- as.formula(paste(outcome,
                            paste(vars, collapse = "+"),
                            sep = "~"))

#Hold all model output in list
rckdri_cox_dcgf<- list()

for(i in 1:n.imp) {
  #Data for imputation i
  tempD<- filter(D, X_mi_m==i)
  
  #Include KDRI LP in a Cox model
  rckdri_cox_dcgf[[i]]<- coxph(kdriform, data = tempD, ties = "efron", x=T)
}

coefs<- data.frame()
for(i in 1:n.imp) {
  coefs<- rbind(coefs, data.frame(imputation=i,
                                  var = vars,
                                  summary(rckdri_cox_dcgf[[i]])$coefficients[,c("coef", "se(coef)")]))
}

#Pool coefficients
pool_coef<- coefs %>%
  group_by(var) %>%
  mutate(pool_coef = mean(coef),
         var_within = sum(se.coef.^2)/n.imp,
         var_between = sum((coef - pool_coef)^2) / (n.imp-1),
         var_total = var_within + var_between + var_between/n.imp,
         pool_se = sqrt(var_total),
         lower = pool_coef - qnorm(0.975)*pool_se,
         upper = pool_coef + qnorm(0.975)*pool_se) %>%
  dplyr::select(-contains("var_")) %>%
  mutate_if(is.numeric, round, 3) %>%
  mutate(ci = paste0("[",lower,", ",upper,"]"),
         out = paste0(pool_coef, "<br>", ci)) %>%
distinct(var, pool_coef, pool_se, ci, out)

write.csv(pool_coef, "../output/recalibration/by_imputation/cox_dcgf_poolcoef.csv", row.names = F)


#Linear predictor
D$LP<- as.matrix(D[,vars]) %*% as.matrix(pool_coef[,"pool_coef"])

#Distribution of LP
ggplot(D, aes(x=LP-uskdri_lp)) +
  geom_histogram() +
  facet_wrap(~X_mi_m)
#Some cases where difference between LP and KDRI is very different
#Does this matter?
#Does this have any affect of why calibration is so different in imputed 7

filter(D, LP-uskdri_lp<=-0.5) %>%
  select(vars, X_mi_m, uskdri_lp, LP) %>%
  View()

#DEFINE DATA FRAMES AND LISTS
rckdrival_cox_dcgf<- list()

#C index
cox_dcgf_cstat<- data.frame()

B<- 500 #Bootstrap samples
cox_dcgf_cboot<- data.frame(matrix(nrow=B)) #Data frame to store C-index from each bootstrap sample

#Store results from score function
cox_dcgf_score<- list()

#Store AUC and Brier scores
cox_dcgf_auc<- data.frame()
cox_dcgf_brier<- data.frame()

#Store Aalen-Johansen and predicted risks
cox_dcgf_risks<- data.frame()



#START LOOP
for (i in 1:n.imp) {
  #Data for imputation i
  tempD<- filter(D, X_mi_m==i)
  
  #Include KDRI LP in a Cox model
  rckdrival_cox_dcgf[[i]]<- coxph(Surv(gsurv_yr, gcens) ~ LP, data = tempD, ties = "efron",
                                  x=T) #x=T is for riskRegression::Score
  
  print("models ran")
  
  #C-index using IPCW because of Wolbers 2014
  #Apparent C index
  cox_dcgf_tempout<- pec::cindex(rckdrival_cox_dcgf[[i]],
                                 formula = Surv(gsurv_yr, gcens) ~ 1,
                                 data = tempD,
                                 splitMethod = "BootCv",
                                 B = B,
                                 keep.matrix = T,
                                 verbose=F)
  print("Store cindex")
  
  cox_dcgf_cstat<- rbind(cox_dcgf_cstat, data.frame(imputation=i,
                                                    appCstat=cox_dcgf_tempout$AppCindex$coxph,
                                                    bootCstat=cox_dcgf_tempout$BootCvCindex$coxph))
  
  print("cindex stored")
  
  cox_dcgf_cboot<- cbind(cox_dcgf_cboot, cox_dcgf_tempout$BootstrapCrossValCindexMat[[1]])
  # #Bootstrap to get optimiasm corrected
  # for (j in 1:B) {
  #   bootdata<- tempD[sample(nrow(tempD), replace = TRUE),]
  #   cox_dcgf_tempboot<- pec::cindex(rckdrival_cox_dcgf[[i]],
  #                                   formula = Surv(gsurv_yr, gcens) ~ 1,
  #                                   data = bootdata)
  #   cox_dcgf_cboot[j,i]<- cox_dcgf_tempboot$AppCindex$coxph
  #   
  #   if(j %in% seq(10,B, by=10)) print(paste0("Imputation ",i,", bootstrap ",j))
  # }
  # cox_dcgf_out[i,"appCstat_sd"]<- sd(cox_dcgf_tempout$BootstrapCrossValCindexMat)
  # cox_dcgf_out[i,"appCstat_lower"]<- quantile(cox_dcgf_tempout$BootstrapCrossValCindexMat, probs=0.25)
  # cox_dcgf_out[i,"appCstat_upper"]<- quantile(cox_dcgf_tempout$BootstrapCrossValCindexMat, probs=0.75)
  # 
  
  print(paste0("Time-dep for imputation", i))
  
  #Time-dependent AUC and Brier score at 1-year and 5-years following transplantation
  cox_dcgf_score[[i]] <- riskRegression::Score(
    list("kdri_val" = rckdrival_cox_dcgf[[i]]),
    formula = Surv(gsurv_yr, gcens) ~ 1, 
    cens.model = "km", 
    data = tempD, 
    conf.int = TRUE, 
    times = c(1,5),
    metrics = c("auc", "brier"),
    summary = c("risks"),
    plots = "calibration",
    split.method = "bootcv",
    B = B)
  
  cox_dcgf_auc<- rbind(cox_dcgf_auc, data.frame(imputation=i,
                                                filter(cox_dcgf_score[[i]]$AUC$score, model=="kdri_val")[,c("times", "AUC", "se")]))
  
  cox_dcgf_brier<- rbind(cox_dcgf_brier, data.frame(imputation=i,
                                                    filter(cox_dcgf_score[[i]]$Brier$score, model=="kdri_val")[,c("times", "Brier", "se")]))

  
  
  #Predicted risks and pseudo values
  cox_dcgf_risks<- rbind(cox_dcgf_risks, data.frame(imputation=i, cox_dcgf_score[[i]]$Calibration$plotframe))
  
  print(paste0("imputation ",i, " done"))
  
}

saveRDS(rckdri_cox_dcgf, file = "../output/rc_cox_dcgf_model.RData")
saveRDS(cox_dcgf_score, file = "../output/rc_cox_dcgf_score.RData")

write.csv(cox_dcgf_auc, file="../output/recalibration/by_imputation/rc_cox_dcgf_auc.csv", row.names = F)
write.csv(cox_dcgf_brier, file="../output/recalibration/by_imputation/rc_cox_dcgf_brier.csv", row.names = F)
write.csv(cox_dcgf_cstat, file="../output/recalibration/by_imputation/rc_cox_dcgf_cstat.csv", row.names = F)
write.csv(cox_dcgf_cboot, file="../output/recalibration/by_imputation/rc_cox_dcgf_cboot.csv", row.names = F)
write.csv(cox_dcgf_risks, file="../output/rc_cox_dcgf_risks.csv", row.names = F)



apply(cox_dcgf_cboot[,-1], 2, quantile, probs=0.75)






#Calibration
n.imp<- 10
cox_dcgf_cal<- data.frame(matrix(nrow = n.imp, ncol = 2))
cox_dcgf_calboot<- data.frame()
B<- 500

rc_cox_dcgf_score<- readRDS("../output/rc_cox_dcgf_score.RData")

timehor<- c(1,5)

for(j in 1:length(timehor)) {
  pdf(file = paste0("../output/recalibration/by_imputation/rc_cox_dcgf_cal_",timehor[j],"yr.pdf"))
  
  for(i in 1:n.imp) {
    #Extract pseudo observations from score function
    pseudos<- filter(rc_cox_dcgf_score[[i]]$Calibration$plotframe, times == timehor[j])
    pseudos_sorted<- arrange(pseudos, risk) #Sort by risk
    
    #Smooth pseudo values using weighted local regression (LOESS)
    loess_pseudo<- predict(loess(pseudovalue ~ risk, data = pseudos_sorted,
                                 degree = 1, #Fit polynomial of degree 1 (linear) between groups
                                 span = 0.2 #Proportion of closest points to use in fit (span*number of obs)
    ),
    se = T)
    
    #Calibration slope
    #Clog-log risk estimates
    pseudos_sorted<- mutate(pseudos_sorted, clog_risks = log(-log(1 - risk)))
    
    #Fit model to estimate calibration slope
    cal_slope<- geepack::geese(pseudovalue ~ offset(clog_risks) + clog_risks,
                               data = pseudos_sorted, id = ID,
                               scale.fix = TRUE,
                               family = gaussian,
                               mean.link = "cloglog",
                               corstr = "independence",
                               jack = TRUE)
    cox_dcgf_cal[i,j]<- summary(cal_slope)$mean["clog_risks","estimate"]
    
    #Bootstrap to get optimism-corrected calibration slope
    for(k in 1:B) {
      bootdata<- pseudos_sorted[sample(nrow(pseudos_sorted), replace = TRUE),]
      tempcal_slope<- geepack::geese(pseudovalue ~ offset(clog_risks) + clog_risks,
                                 data = bootdata, id = ID,
                                 scale.fix = TRUE,
                                 family = gaussian,
                                 mean.link = "cloglog",
                                 corstr = "independence",
                                 jack = TRUE)
      cox_dcgf_calboot<- rbind(cox_dcgf_calboot, data.frame(imputation=i,
                                                            times=timehor[j],
                                                            calslope=summary(tempcal_slope)$mean["clog_risks","estimate"],
                                                            cal_se= summary(tempcal_slope)$mean["clog_risks","san.se"]))
    }
    
    
    #PLOT
    #Setting up
    axislim<- ceiling(max(loess_pseudo$fit + qt(p = 0.975, df = loess_pseudo$df) * loess_pseudo$se.fit)/0.1)*0.1
    spike_bounds <- c(-0.05, 0)
    bin_breaks <- seq(0, axislim, length.out = 100 + 1)
    freqs <- table(cut(pseudos_sorted$risk, breaks = bin_breaks))
    bins <- bin_breaks[-1]
    freqs_valid <- freqs[freqs > 0]
    freqs_rescaled <- spike_bounds[1] + (spike_bounds[2] - spike_bounds[1]) * 
      (freqs_valid - min(freqs_valid)) / (max(freqs_valid) - min(freqs_valid))
    
    #Start with a blank plot
    plot(x = pseudos_sorted$risk, y = pseudos_sorted$pseudovalue,
         xlim = c(0, axislim), ylim = c(spike_bounds[1], axislim),
         yaxt = "n", xaxs="i",
         xlab = "Estimated risks", ylab = "Observed outcome proprtion",
         frame.plot = FALSE, type = "n")
    axis(side = 2, at = seq(0, axislim, by=0.1))
    #Add confidence intervals
    polygon(x = c(pseudos_sorted$risk, rev(pseudos_sorted$risk)),
            y = c(pmax(loess_pseudo$fit - qt(p = 0.975, df = loess_pseudo$df) * loess_pseudo$se.fit, 0), #Stop confidence interval dipping below 0
                  rev(loess_pseudo$fit + qt(p = 0.975, df = loess_pseudo$df) * loess_pseudo$se.fit)),
            col = "lightgrey")
    #Add diagonal line for reference of perfect calibration
    abline(a=0, b=1, col="red", lty=2)
    #Add line of calibration
    lines(x = pseudos_sorted$risk, y = loess_pseudo$fit,
          col = "black", lwd = 2)
    #Add histogram of predicted risk underneath x axis
    segments(x0 = bins[freqs > 0], x1 = bins[freqs > 0],
             y0 = spike_bounds[1], y1 = freqs_rescaled)
    
    title(paste0("Imputation ",i,", ",timehor[j],"-year death-censored graft failure"))
    
    print(i)
  }
  dev.off()
}

write.csv(cox_dcgf_calboot, "../output/recalibration/by_imputation/cox_dcgf_calboot.csv", row.names = F)
write.csv(cox_dcgf_cal, "../output/recalibration/by_imputation/cox_dcgf_cal.csv", row.names = F)


cox_dcgf_calboot %>%
  group_by(times) %>%
  summarise(mean=mean(calslope),
            lower=quantile(calslope, 0.25),
            upper=quantile(calslope, 0.75),
            se=(upper-lower)/(2*qnorm(0.975)))
  
filter(cox_dcgf_calboot, times==1) %>%
  ggplot(aes(x=calslope)) +
  geom_histogram() +
  facet_wrap(~imputation)


