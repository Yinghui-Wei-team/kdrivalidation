################################################
###             KDRI recalibration           ###
###           Composite outcome model        ###
################################################

#Packages
library(tidyverse)
library(prodlim)
library(survival)
library(riskRegression)
library(pec)
library(geepack)




# Data --------------------------------------------------------------------
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



# Initialise --------------------------------------------------------------
#Recalibrate the KDRI using the cause-specific Cox model
#Number of imputations
n.imp<- 15

#Bootstrap samples
B<- 100

#Event time horizons
times<- c(1,5)

#Model formula
outcome<- "Surv(gsurv_yr, gcens)"
vars<- c("rcdage18", "rcdage", "rcdage50", "rcdheight", "rcdweight",
         "dethnicb_i", "dpast_hypertension_i", "dpast_diabetes_i", "cva_i",
         "rcdcreat", "rcdcreat15", "dhcv_i", "dtype_i")

kdriform<- as.formula(paste(outcome,
                            paste(vars, collapse = "+"),
                            sep = "~"))

#Hold all model output in list
rckdri_dcgf<- list()

#Coefficients from regression model
coefs_dcgf<- data.frame()

#Store results from score function
rc_dcgf_score_app<- list()

#Store predicted risks
rc_dcgf_risks<- data.frame()

#Store Numeric performance measures
rc_dcgf_auc<- data.frame()
rc_dcgf_brier<- data.frame()
rc_dcgf_cal<- data.frame()





# START LOOP --------------------------------------------------------------

for (k in 1:length(times)) {
  
  #Save calibration plots as pdf
  pdf(file = paste0("../output/recalibration/by_imputation/dcgf_cal_",times[k],"yr.pdf"))
  
  for (i in 1:n.imp) {
    #Data for imputation i
    tempD<- filter(D, X_mi_m==i)
    
    ## Fit model ----------------------------------------------------------
    #Include KDRI variables in a cause-specific Cox model
    rckdri_dcgf[[i]]<- coxph(kdriform, data = tempD, ties = "efron",
                             x=T) #x=T is for riskRegression::Score
    
    #Store coefficients -> to be pooled afterwards
    coefs_dcgf<- rbind(coefs_dcgf, data.frame(imputation=i,
                                              var = vars,
                                              summary(rckdri_dcgf[[i]])$coefficients[,c("coef", "se(coef)")]))
    

    
    ## Apparent performance -----------------------------------------------
    #Discrimination
    #Time-dependent AUC and Brier score at 1-year and 5-years following transplantation
    rc_dcgf_score_app[[i]] <- riskRegression::Score(
      list("kdri_rc" = rckdri_dcgf[[i]]),
      formula = Surv(gsurv_yr, gcens) ~ 1,
      cens.model = "km", 
      data = tempD, 
      conf.int = TRUE, 
      times = times[k],
      metrics = c("auc", "brier"),
      summary = c("risks"),
      plots = "calibration")
    
    
    
    #Predicted risk from model and pseudo values
    rc_dcgf_risks<- rbind(rc_dcgf_risks, data.frame(imputation=i, 
                                                    rc_dcgf_score_app[[i]]$Calibration$plotframe))
    
    
    
    #Calibration
    #Extract pseudo observations from score function
    pseudos<- rc_dcgf_score_app[[i]]$Calibration$plotframe
    pseudos_sorted<- arrange(pseudos, risk) #Sort by risk
    
    #Smooth pseudo values using weighted local regression (LOESS)
    loess_pseudo<- predict(loess(pseudovalue ~ risk, data = pseudos_sorted,
                                 degree = 1, #Fit polynomial of degree 1 (linear) between groups
                                 span = 0.2 #Proportion of closest points to use in fit (span*number of obs)
    ),
    se = T)
    
    print("loess smoothing done")
    
    #Calibration plot
    #NOTE: The calibration plot represents the apparent calibration
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
         main = paste0("Imputation ",i,", ",times[k]," year death-censored \n graft failure"))
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
    
    
    
    #Calibration slope
    #Clog-log risk estimates
    pseudos_sorted<- mutate(pseudos_sorted, clog_risks = log(-log(1 - risk)))
    
    #Fit model to estimate calibration slope
    #NOTE: This should equal 1
    cal_slope_app<- geepack::geese(pseudovalue ~ offset(clog_risks) + clog_risks,
                                   data = pseudos_sorted, id = ID,
                                   scale.fix = TRUE,
                                   family = gaussian,
                                   mean.link = "cloglog",
                                   corstr = "independence"#,
                                   #jack = TRUE
                                   )
    
    
    
    
    
    print(paste0("Apparent performance for imputation ", i, " in time ", times[k], " complete!"))
    
    
    ## Optimism-adjusted performance --------------------------------------
    for (j in 1:B) {
      #Bootstrapped data
      boot<- sample(nrow(tempD), replace = TRUE)
      bootdata<- tempD[boot, ]
      
      #Fit  model to bootstrapped data
      rckdri_dcgf_boot<- coxph(kdriform, data = tempD, ties = "efron",
                               x=T) #x=T is for riskRegression::Score
      
      
      #Discrimination
      #Performance of the bootstrapped model (AUC_{i,boot}, brier_{i,boot})
      rc_dcgf_score_boot<- riskRegression::Score(
        list("kdri_rc_boot" = rckdri_dcgf_boot),
        formula = Surv(gsurv_yr, gcens) ~ 1, 
        cens.model = "km", 
        data = tempD[sample(nrow(tempD), replace = TRUE),], 
        conf.int = FALSE, 
        times = times[k],
        metrics = c("auc", "brier"),
        plots = "calibration")
      
      #Performance of bootstrapped model in original (full) data (AUC_{i,orig}, brier_{i,orig})
      rc_dcgf_score_orig<- riskRegression::Score(
        list("kdri_rc_boot" = rckdri_dcgf_boot),
        formula = Surv(gsurv_yr, gcens) ~ 1, 
        cens.model = "km", 
        data = tempD, 
        conf.int = FALSE, 
        times = times[k],
        metrics = c("auc", "brier"),
        plots = "calibration")
      
      #Store results
      rc_dcgf_auc<- rbind(rc_dcgf_auc, data.frame(imputation=i,
                                                boot=j,
                                                time = times[k],
                                                auc_app = filter(rc_dcgf_score_app[[i]]$AUC$score, model=="kdri_rc")[,c("AUC")],
                                                auc_boot = filter(rc_dcgf_score_boot$AUC$score, model=="kdri_rc_boot")[,c("AUC")],
                                                auc_orig = filter(rc_dcgf_score_orig$AUC$score, model=="kdri_rc_boot")[,c("AUC")]
      )
      )
      
      rc_dcgf_brier<- rbind(rc_dcgf_brier, data.frame(imputation=i,
                                                    boot=j,
                                                    time = times[k],
                                                    brier_app = filter(rc_dcgf_score_app[[i]]$Brier$score, model=="kdri_rc")[,c("Brier")],
                                                    brier_boot = filter(rc_dcgf_score_boot$Brier$score, model=="kdri_rc_boot")[,c("Brier")],
                                                    brier_orig = filter(rc_dcgf_score_orig$Brier$score, model=="kdri_rc_boot")[,c("Brier")]
      )
      )
      
      
      #Calibration
      #Performance of the bootstrapped model (CalSlope_{i,boot})
      #Extract pseudo observations from score function
      pseudos_boot<- arrange(rc_dcgf_score_boot$Calibration$plotframe, risk)
      
      
      #Calibration slope
      #Clog-log risk estimates
      pseudos_boot<- mutate(pseudos_boot, clog_risks = log(-log(1 - risk)))
      
      #Fit model to estimate calibration slope
      cal_slope_boot<- geepack::geese(pseudovalue ~ offset(clog_risks) + clog_risks,
                                      data = pseudos_boot, id = ID,
                                      scale.fix = TRUE,
                                      family = gaussian,
                                      mean.link = "cloglog",
                                      corstr = "independence"#,
                                      #jack = TRUE
                                      )
      
      
      #Performance of the bootstrapped model in original (full) data (CalSlope_{i,full})
      #Extract pseudo observations from score function
      pseudos_orig<- arrange(rc_dcgf_score_orig$Calibration$plotframe, risk)
      
      #Clog-log risk estimates
      pseudos_orig<- mutate(pseudos_orig, clog_risks = log(-log(1 - risk)))
      
      #Fit model to estimate calibration slope
      cal_slope_orig<- geepack::geese(pseudovalue ~ offset(clog_risks) + clog_risks,
                                      data = pseudos_orig, id = ID,
                                      scale.fix = TRUE,
                                      family = gaussian,
                                      mean.link = "cloglog",
                                      corstr = "independence"#,
                                      #jack = TRUE
                                      )
      
      #Store results
      rc_dcgf_cal<- rbind(rc_dcgf_cal, data.frame(imputation=i,
                                                boot=j,
                                                time = times[k],
                                                cal_app = 1 + summary(cal_slope_app)$mean["clog_risks","estimate"],
                                                cal_boot = 1 + summary(cal_slope_boot)$mean["clog_risks","estimate"],
                                                cal_orig = 1 + summary(cal_slope_orig)$mean["clog_risks","estimate"]
      )
      )
      
    } #END OF BOOTSTRAP
    
    #print(paste0("Optimism-corrected performance for imputation ", i, " in time ", times[k], " complete!"))
    
    print(paste0("End imputation ", i, " for time horizon ", times[k], ": ", Sys.time()))
  }
  
  #Start new plot for new time horizon
  dev.off()
  
}


# Write out results -------------------------------------------------------

#Model
write.csv(coefs_dcgf, "../output/recalibration/by_imputation/coefs_dcgf.csv", row.names = F)
saveRDS(rckdri_dcgf, file = "../output/rc_dcgf_model.RData")

#Predicted risks
write.csv(rc_dcgf_risks, "../output/recalibration/by_imputation/rc_dcgf_risks.csv", row.names = F)

#Performance measures
colnames(rc_dcgf_auc)<- c("imputation", "boot", "time", "auc_app", "auc_boot", "auc_orig")
colnames(rc_dcgf_brier)<- c("imputation", "boot", "time", "brier_app", "brier_boot", "brier_orig")
colnames(rc_dcgf_cal)<- c("imputation", "boot", "time", "cal_app", "cal_boot", "cal_orig")

write.csv(rc_dcgf_auc, file="../output/recalibration/by_imputation/rc_dcgf_auc.csv", row.names = F)
write.csv(rc_dcgf_brier, file="../output/recalibration/by_imputation/rc_dcgf_brier.csv", row.names = F)
write.csv(rc_dcgf_cal, file="../output/recalibration/by_imputation/rc_dcgf_cal.csv", row.names = F)












