################################################
###             KDRI recalibration           ###
###            Competing risk model          ###
################################################

#Packages
library(tidyverse)
library(prodlim)
library(survival)
library(riskRegression)
library(pec)



#Data
data<- read.csv("../../NHSBT Data/Data/kdri_mi_15.csv", na.strings = "")

D<- data %>%
  filter(X_mi_m > 0 & X_mi_m <=15) %>%
  mutate(gsurv = case_when(gsurv==0 ~ 0.5,
                           TRUE ~ as.numeric(gsurv)),
         psurv = case_when(psurv==0 ~ 0.5,
                           TRUE~ as.numeric(psurv)),
         gsurv_yr=gsurv/365.5,
         psurv_yr=psurv/365.5) %>%
  select(arecip_id, gsurv_yr, gcens, psurv_yr, pcens, uskdri_lp, X_mi_m,
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


D_cr<- D %>%
  mutate(etime = case_when(gcens==1 ~ gsurv_yr,
                           TRUE ~ psurv_yr),
         etype = case_when(gcens==0 ~ 2*pcens,
                           TRUE ~ 1),
         etype_fct = factor(etype))


#Recalibrate the KDRI using the cause-specific Cox model
#Number of imputations
n.imp<- 15

#Model formula
outcome<- "Hist(etime, etype_fct)"
vars<- c("rcdage18", "rcdage", "rcdage50", "rcdheight", "rcdweight",
         "dethnicb_i", "dpast_hypertension_i", "dpast_diabetes_i", "cva_i",
         "rcdcreat", "rcdcreat15", "dhcv_i", "dtype_i")

kdriform<- as.formula(paste(outcome,
                            paste(vars, collapse = "+"),
                            sep = "~"))

#Update by developing with a competing risks model
#Hold all model output in list
rckdri_csc<- list()

B<- 500 #Bootstrap samples

#Store results from score function
csc_score_app<- list()
csc_score_boot<- list()
csc_score_orig<- list()


#Store AUC and Brier scores
csc_auc<- data.frame()
csc_brier<- data.frame()

#Store Aalen-Johansen and predicted risks
csc_risks<- data.frame()



#START LOOP
for (i in 1:n.imp) {
  #Data for imputation i
  tempD<- filter(D_cr, X_mi_m==i)
  
  #Include KDRI variables in a cause-specific Cox model
  rckdri_csc[[i]]<- CSC(kdriform, data = tempD,
                        cause = 1)
  
  print("model ran")
  
  
  #Time-dependent AUC and Brier score at 1-year and 5-years following transplantation
  csc_score_app[[i]] <- riskRegression::Score(
    list("kdri_rc" = rckdri_csc[[i]]),
    formula = Hist(etime, etype_fct) ~ 1, 
    cause = 1,
    cens.model = "km", 
    data = tempD, 
    conf.int = TRUE, 
    times = c(1,5),
    metrics = c("auc", "brier"),
    summary = c("risks"),
    plots = "calibration")
  
  print("apparent perf done")
  
  #Harrell's optimism corrected performance measures
  for (j in 1:B) {
    #Construct model for B bootstrap samples
    rckdri_csc_boot<- CSC(kdriform, data = tempD[sample(nrow(tempD), replace = TRUE),],
                          cause = 1)
    
    #Performance of the bootstrapped model (AUC_{i,boot}, brier_{i,boot})
    csc_score_boot<- riskRegression::Score(
      list("kdri_rc" = rckdri_csc_boot),
      formula = Hist(etime, etype_fct) ~ 1, 
      cause = 1,
      cens.model = "km", 
      data = tempD[sample(nrow(tempD), replace = TRUE),], 
      conf.int = FALSE, 
      times = c(1,5),
      metrics = c("auc", "brier"))
    
    #Performance of bootstrapped model in original (full) data (AUC_{i,orig}, brier_{i,orig})
    csc_score_orig<- riskRegression::Score(
      list("kdri_rc" = rckdri_csc_boot),
      formula = Hist(etime, etype_fct) ~ 1, 
      cause = 1,
      cens.model = "km", 
      data = tempD, 
      conf.int = FALSE, 
      times = c(1,5),
      metrics = c("auc", "brier"))
    
    #Store results
    csc_auc<- rbind(csc_auc, data.frame(imputation=i,
                                        boot=j,
                                        time = c(1,5),
                                        auc_app = filter(csc_score_app[[i]]$AUC$score, model=="kdri_rc")[,c("AUC")],
                                        auc_boot = filter(csc_score_boot$AUC$score, model=="kdri_rc")[,c("AUC")],
                                        auc_orig = filter(csc_score_orig$AUC$score, model=="kdri_rc")[,c("AUC")]
                                        )
                    )
    
    # csc_brier<- rbind(csc_brier, data.frame(imputation=i,
    #                                             filter(csc_score[[i]]$Brier$score, model=="kdri_rc")[,c("Brier")]))
    
    csc_brier<- rbind(csc_brier, data.frame(imputation=i,
                                            boot=j,
                                            time = c(1,5),
                                            brier_app = filter(csc_score_app[[i]]$Brier$score, model=="kdri_rc")[,c("Brier")],
                                            brier_boot = filter(csc_score_boot$Brier$score, model=="kdri_rc")[,c("Brier")],
                                            brier_orig = filter(csc_score_orig$Brier$score, model=="kdri_rc")[,c("Brier")]
                                            )
                      )
    
    
  }
  

  
  #Predicted risk from model and pseudo values
  csc_risks<- rbind(csc_risks, data.frame(imputation=i, csc_score_app[[i]]$Calibration$plotframe))
  
  print(paste0("imputation ",i, " done"))
  
  
}

saveRDS(rckdri_csc, file = "../output/rc_csc_model.RData")
saveRDS(csc_score_app, file = "../output/rc_csc_score.RData")

colnames(csc_auc)<- c("imputation", "boot", "time", "auc_app", "auc_boot", "auc_orig")
colnames(csc_brier)<- c("imputation", "boot", "time", "brier_app", "brier_boot", "brier_orig")

write.csv(csc_auc, file="../output/recalibration/by_imputation/rc_csc_auc.csv", row.names = F)
write.csv(csc_brier, file="../output/recalibration/by_imputation/rc_csc_brier.csv", row.names = F)
write.csv(csc_risks, file="../output/rc_csc_risks.csv", row.names = F)


#Calibration
csc_cal<- data.frame(matrix(nrow = n.imp, ncol = 2))
csc_calboot<- data.frame()

rc_csc_score<- readRDS("../output/rc_csc_score.RData")
timehor<- c(1,5)

print(paste0("Start time:", Sys.time()))

for(j in 1:length(timehor)) {
  pdf(file = paste0("../output/recalibration/by_imputation/csc_cal_",timehor[j],"yr.pdf"))
  
  for(i in 1:n.imp) {
    #Extract pseudo observations from score function
    pseudos<- filter(rc_csc_score[[i]]$Calibration$plotframe, times == timehor[j])
    pseudos_sorted<- arrange(pseudos, risk) #Sort by risk
    
    #Smooth pseudo values using weighted local regression (LOESS)
    loess_pseudo<- predict(loess(pseudovalue ~ risk, data = pseudos_sorted,
                                 degree = 1, #Fit polynomial of degree 1 (linear) between groups
                                 span = 0.2 #Proportion of closest points to use in fit (span*number of obs)
    ),
    se = T)
    
    print("loess smoothing done")
    
    #Calibration slope
    #Clog-log risk estimates
    pseudos_sorted<- mutate(pseudos_sorted, clog_risks = log(-log(1 - risk)))
    
    #Fit model to estimate calibration slope
    cal_slope_app<- geepack::geese(pseudovalue ~ offset(clog_risks) + clog_risks,
                               data = pseudos_sorted, id = ID,
                               scale.fix = TRUE,
                               family = gaussian,
                               mean.link = "cloglog",
                               corstr = "independence",
                               jack = TRUE)
    csc_cal[i,j]<- summary(cal_slope)$mean["clog_risks","estimate"]
    
    print("start bootstrap for cal slope")
    
    #Bootstrap to get optimism-corrected calibration slope
    for(k in 1:B) {
      #Bootstrap samples
      bootdata<- pseudos_sorted[sample(nrow(pseudos_sorted), replace = TRUE),]
      
      #Calibration slope in bootstrap data (cal_)
      cal_slope_boot<- geepack::geese(pseudovalue ~ offset(clog_risks) + clog_risks,
                                     data = bootdata, id = ID,
                                     scale.fix = TRUE,
                                     family = gaussian,
                                     mean.link = "cloglog",
                                     corstr = "independence",
                                     jack = TRUE)
      
      #Calibration slope in 
      cal_slope_boot<- geepack::geese(pseudovalue ~ offset(clog_risks) + clog_risks,
                                      data = bootdata, id = ID,
                                      scale.fix = TRUE,
                                      family = gaussian,
                                      mean.link = "cloglog",
                                      corstr = "independence",
                                      jack = TRUE)
      
      #Store results
      csc_calboot<- rbind(csc_calboot, data.frame(imputation=i,
                                                  time=timehor[j],
                                                  calslope_app = summary(cal_slope_app)$mean["clog_risks","estimate"],
                                                  calslope_boot = summary(cal_slope_boot)$mean["clog_risks","estimate"],
                                                  calslope_orig = summary(cal_slope_orig)$mean["clog_risks","estimate"]
                                                  )
                          )
    }
    
    
    #PLOT with loess smoothing
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
         frame.plot = FALSE, type = "n",
         main = paste0("Imputation ",i,", ",timehor[j]," year graft failure \n (death as competing event)"))
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
    
    print(i)
    
    print(paste0("End imputation ", i, " for time horizon ", timehor[j], ":", Sys.time()))
  }
  dev.off()
}

write.csv(csc_calboot, "../output/recalibration/by_imputation/csc_calboot.csv", row.names = F)
write.csv(csc_cal, "../output/recalibration/by_imputation/csc_cal.csv", row.names = F)


csc_calboot %>%
  group_by(times) %>%
  summarise(mean=mean(calslope),
            lower=quantile(calslope, 0.25),
            upper=quantile(calslope, 0.75),
            se=(upper-lower)/(2*qnorm(0.975)))

filter(csc_calboot, times==1) %>%
  ggplot(aes(x=calslope)) +
  geom_histogram() +
  facet_wrap(~imputation)
