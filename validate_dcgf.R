################################################
###              KDRI validation             ###
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
  select(arecip_id, gsurv_yr, gcens, uskdri_lp, X_mi_m)

# keep arecip_id tsurv tcens gsurv gcens psurv pcens dage* dheight* dweight* dethnic* dpast_hypertension dpast_diabetes dcreat* dhcv dtype
# uskdri_lp = 0.0128*(dage-40) + dage50_i*0.0107*(dage-50) - 0.0464*(dheight-170)/10 - dweight80_i*0.0199*(dweight-80)/5 + 0.179*dethnicb_i + 0.126*dpast_hypertension + 0.13*dpast_diabetes + 0.0881*cva_i + 0.22*(dcreatmgdl-1) - dcreat15_i*0.209*(dcreatmgdl-1.5) + 0.24*dhcv + 0.133*dtype))



#Validate using the Cox model
#Number of imputations
n.imp<- 10

#Hold all model output in list
kdri_cox_dcgf<- list()

#C index
cox_dcgf_out<- data.frame(matrix(nrow=n.imp, ncol = 9))
colnames(cox_dcgf_out)<- c("imputation", "Cstat", "Cstat_sd", "Cstat_lower", "Cstat_upper", "Cslope", "Cslope_se", "Cslope_lower", "Cslope_upper")
B<- 100 #Bootstrap samples
cox_dcgf_cboot<- c() #Vector to store C-index from each bootstrap sample

#Store results from score function
cox_dcgf_score<- list()

#Store AUC and Brier scores
cox_dcgf_auc_brier<- data.frame()

#Store Aalen-Johansen and predicted risks
cox_dcgf_risks<- data.frame()



#START LOOP
for (i in 1:n.imp) {
  #Data for imputation i
  tempD<- filter(D, X_mi_m==i)
  
  #Include KDRI LP in a Cox model
  kdri_cox_dcgf[[i]]<- coxph(Surv(gsurv_yr, gcens) ~ uskdri_lp, data = tempD, ties = "efron",
                             x=T) #x=T is for riskRegression::Score
  
  
  #C-index using IPCW because of Wolbers 2014
  #Estimate
  cox_dcgf_out[i,1]<- i
  cox_dcgf_out[i,2]<- pec::cindex(kdri_cox_dcgf[[i]],
                                  formula = Surv(gsurv_yr, gcens) ~ 1,
                                  data = tempD)$AppCindex$coxph
  #Bootstrap to get SE of C-index
  for (j in 1:B) {
    bootdata<- tempD[sample(nrow(tempD), replace = TRUE),]
    cox_dcgf_cboot[j]<- pec::cindex(kdri_cox_dcgf[[i]],
                                    formula = Surv(gsurv_yr, gcens) ~ 1,
                                    data = bootdata)$AppCindex$coxph
    if(j %in% c(10,20,30,40,50,60,70,80,90)) print(paste0("Imputation ",i,", bootstrap ",j))
  }
  cox_dcgf_out[i,3]<- sd(cox_dcgf_cboot)
  cox_dcgf_out[i,4]<- quantile(cox_dcgf_cboot, probs=0.25)
  cox_dcgf_out[i,5]<- quantile(cox_dcgf_cboot, probs=0.75)
  
  
  #Calibration slope
  cox_dcgf_out[i, 6]<- kdri_cox_dcgf[[i]]$coefficients
  cox_dcgf_out[i, 7]<- kdri_cox_dcgf[[i]]$var
  cox_dcgf_out[i, 8]<- kdri_cox_dcgf[[i]]$coefficients - (qnorm(0.975)*kdri_cox_dcgf[[i]]$var)
  cox_dcgf_out[i, 9]<- kdri_cox_dcgf[[i]]$coefficients + (qnorm(0.975)*kdri_cox_dcgf[[i]]$var)
  
  
  #Time-dependent AUC and Brier score at 1-year and 5-years following transplantation
  cox_dcgf_score[[i]] <- riskRegression::Score(
    list("kdri_val" = kdri_cox_dcgf[[i]]),
    formula = Surv(gsurv_yr, gcens) ~ 1, 
    cens.model = "km", 
    data = tempD, 
    conf.int = TRUE, 
    times = c(1,5),
    metrics = c("auc", "brier"),
    summary = c("risks"),
    plots = "calibration"
  )
  cox_dcgf_auc_brier<- rbind(cox_dcgf_auc_brier, data.frame(imputation=i, 
                                                            cox_dcgf_score[[i]]$AUC$score[,c("times", "AUC", "se", "lower", "upper")], 
                                                            filter(cox_dcgf_score[[i]]$Brier$score, model=="kdri_val")[,c("times", "Brier", "se", "se.conservative", "lower", "upper")]))
  
  
  #Predicted risks and pseudo values
  cox_dcgf_risks<- rbind(cox_dcgf_risks, data.frame(imputation=i, cox_dcgf_score[[i]]$Calibration$plotframe))
  
  
}

saveRDS(kdri_cox_dcgf, file = "../output/cox_dcgf_model.RData")
saveRDS(cox_dcgf_score, file = "../output/cox_dcgf_score.RData")



write.csv(cox_dcgf_auc_brier, file="../output/by_imputation/cox_dcgf_auc_brier.csv", row.names = F)
# write.csv(cox_dcgf_out, file="output/by_imputation/cox_dcgf_cind_slope.csv", row.names = F)
write.csv(cox_dcgf_risks, file="../output/cox_dcgf_risks.csv", row.names = F)



#Calibration plot LOESS smoothing
cox_dcgf_score<- readRDS("../output/cox_dcgf_score.RData")
cox_dcgf_cal<- data.frame()
timehor<- c(1,5)
for(j in 1:length(timehor)) {
  pdf(file = paste0("output/by_imputation/cox_dcgf_cal_",timehor[j],"yr.pdf"))
  
  for(i in 1:n.imp) {
    #Extract pseudo observations from score function
    pseudos<- filter(cox_dcgf_score[[i]]$Calibration$plotframe, times == timehor[j])
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
    cox_dcgf_cal<- rbind(cox_dcgf_cal, data.frame(imputation=i,
                                                  times=timehor[j],
                                                  calslope=summary(cal_slope)$mean["clog_risks","estimate"],
                                                  cal_se= summary(cal_slope)$mean["clog_risks","san.se"]))

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

write.csv(cox_dcgf_cal, "output/by_imputation/cox_dcgf_cal.csv", row.names = F)





#Calibration plots
for(i in c(1,5)) {
  pdf(file = paste0("output/by_imputation/cox_dcgf_cal_",i,"yr.pdf"))
  
  axislim<- ifelse(i==1, 0.2, 0.4)

  for(j in 1:n.imp) {
    plotCalibration(
      x = cox_dcgf_score[[j]],
      times = i,
      bandwidth = 0.03,
      cens.method = "jackknife",
      pseudo = T,
      round = FALSE, 
      xlim = c(0, axislim), ylim = c(0, axislim), 
      rug = TRUE, 
      xlab = "Estimated risks", ylab = "Observed risks",
      auc.in.legend = F, brier.in.legend = F)
    title(paste0("Imputation ",j, ", ", i, "year graft survival (death-censored)"))
  }
dev.off()
}
