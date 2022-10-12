################################################
###              KDRI validation             ###
###            Competing risk model          ###
################################################

#Packages
library(tidyverse)
library(prodlim)
library(survival)
library(riskRegression)
library(pec)



#Data
data<- read.csv("../../NHSBT Data/Data/kdri_mi.csv", na.strings = "")

D<- data %>%
  filter(X_mi_m > 0 & X_mi_m <=10) %>%
  mutate(gsurv = case_when(gsurv==0 ~ 0.5,
                           TRUE ~ as.numeric(gsurv)),
         psurv = case_when(psurv==0 ~ 0.5,
                           TRUE~ as.numeric(psurv)),
         gsurv_yr=gsurv/365.5,
         psurv_yr=psurv/365.5) %>%
  select(arecip_id, gsurv_yr, gcens, psurv_yr, pcens, uskdri_lp, X_mi_m)

D_cr<- D %>%
  mutate(etime = case_when(gcens==1 ~ gsurv_yr,
                           TRUE ~ psurv_yr),
         etype = case_when(gcens==0 ~ 2*pcens,
                           TRUE ~ 1),
         etype_fct = factor(etype))


# D_crprep1<- mstate::crprep(Tstop = c("gsurv_yr","psurv_yr"), status = c("gcens","pcens"),
#                           data = D[which(D$X_mi_m==1),],
#                           trans = 1, cens = 0,
#                           id = "arecip_id",
#                           keep = c("uskdri_lp", "X_mi_m"))



#Validate using a competing risks model
#Number of imputed datasets
n.imp<-10

#Hold all model output in list
kdri_csc<- list()

#C index
csc_out<- data.frame(matrix(nrow=n.imp, ncol = 9))
colnames(csc_out)<- c("imputation", "Cstat", "Cstat_sd", "Cstat_lower", "Cstat_upper", "Cslope", "Cslope_se", "Cslope_lower", "Cslope_upper")
B<- 100 #Bootstrap samples
csc_cboot<- c() #Vector to store C-index from each bootstrap sample

#Store results from score function
csc_score<- list()

#Store AUC and Brier scores
csc_auc_brier<- data.frame()

#Store Aalen-Johansen and predicted risks
csc_risks<- data.frame()



#START LOOP
for (i in 1:n.imp) {
  #Data for imputation i
  tempD<- filter(D_cr, X_mi_m==i)
  
  #Include KDRI LP in a cause-specific Cox model
  kdri_csc[[i]]<- CSC(Hist(etime, etype_fct) ~ uskdri_lp, data = tempD,
                           cause = 1)
  
  
  #C-index using methods of Wolbers 2009
  #Estimate
  csc_out[i,1]<- i
  csc_out[i,2]<- pec::cindex(kdri_csc[[i]],
                             formula = Hist(etime, etype_fct) ~ 1,
                             data = tempD,
                             cause = 1)$AppCindex$CauseSpecificCox

  #Bootstrap to get SE of C-index
  for (j in 1:B) {
    bootdata<- tempD[sample(nrow(tempD), replace = TRUE),]
    csc_cboot[j]<- pec::cindex(kdri_csc[[i]],
                               formula = Hist(etime, etype_fct) ~ 1,
                               data = bootdata,
                               cause = 1)$AppCindex$CauseSpecificCox
    if(j %in% c(10,20,30,40,50,60,70,80,90)) print(paste0("Imputation ",i,", bootstrap ",j))
  }
  csc_out[i,3]<- sd(csc_cboot)
  csc_out[i,4]<- quantile(csc_cboot, probs=0.25)
  csc_out[i,5]<- quantile(csc_cboot, probs=0.75)

  
  #Calibration slope
  csc_out[i, 6]<- kdri_csc[[i]]$models$`Cause 1`$coefficients
  csc_out[i, 7]<- kdri_csc[[i]]$models$`Cause 1`$var
  csc_out[i, 8]<- kdri_csc[[i]]$models$`Cause 1`$coefficients - (qnorm(0.975)*kdri_csc[[i]]$models$`Cause 1`$var)
  csc_out[i, 9]<- kdri_csc[[i]]$models$`Cause 1`$coefficients + (qnorm(0.975)*kdri_csc[[i]]$models$`Cause 1`$var)

  
  #Time-dependent AUC and Brier score at 1-year and 5-years following transplantation
  csc_score[[i]] <- riskRegression::Score(
    list("kdri_val" = kdri_csc[[i]]),
    formula = Hist(etime, etype_fct) ~ 1, 
    cause = 1,
    cens.model = "km", 
    data = tempD, 
    conf.int = TRUE, 
    times = c(1,5),
    metrics = c("auc", "brier"),
    summary = c("risks"),
    plots = "calibration"
  )
  csc_auc_brier<- rbind(csc_auc_brier, data.frame(imputation=i,
                                                  csc_score[[i]]$AUC$score[,c("times", "AUC", "se", "lower", "upper")], 
                                                  filter(csc_score[[i]]$Brier$score, model=="kdri_val")[,c("times", "Brier", "se", "se.conservative", "lower", "upper")]))
  
  
  #Predicted risk from model and pseudo values
  csc_risks<- rbind(csc_risks, data.frame(imputation=i, csc_score[[i]]$Calibration$plotframe))
  
  
}

saveRDS(kdri_csc, file = "../output/csc_model.RData")
saveRDS(csc_score, file = "../output/csc_score.RData")


write.csv(csc_auc_brier, file="output/by_imputation/csc_auc_brier.csv", row.names = F)
write.csv(csc_out, file="output/by_imputation/csc_cind_slope.csv", row.names = F)
write.csv(csc_risks, file="../output/csc_risks.csv", row.names = F)


#Calibration plot LOESS smoothing
csc_score<- readRDS("../output/csc_score.RData")
csc_cal<- data.frame()
timehor<- c(1,5)
for(j in 1:length(timehor)) {
  pdf(file = paste0("output/by_imputation/csc_cal_",timehor[j],"yr.pdf"))
  
  for(i in 1:n.imp) {
    #Extract pseudo observations from score function
    pseudos<- filter(csc_score[[i]]$Calibration$plotframe, times == timehor[j])
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
    cal_slope<- geepack::geese(pseudovalue ~ offset(clog_risks) + clog_risks,
                               data = pseudos_sorted, id = ID,
                               scale.fix = TRUE,
                               family = gaussian,
                               mean.link = "cloglog",
                               corstr = "independence",
                               jack = TRUE)
    csc_cal<- rbind(csc_cal, data.frame(imputation=i,
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
    #Add calibration curve
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

write.csv(csc_cal, "output/by_imputation/csc_cal.csv", row.names = F)


#Number of patients with higher than 10% estimated risk
filter(csc_score[[1]]$Calibration$plotframe, times == 1 & risk>0.1) %>%
  nrow()


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

#Calibration slope
csc_cal<- data.frame()
csc_cal[1,1]<- summary(cal_slope)$mean["clog_risks","estimate"]
csc_cal[1,2]<- summary(cal_slope)$mean["clog_risks","san.se"]



#Calibration plots
for(i in c(1,5)) {
  pdf(file = paste0("output/by_imputation/csc_cal_",i,"yr.pdf"))
  
  axislim<- ifelse(i==1, 0.2, 0.4)
  
  for(j in 1:n.imp) {
    plotCalibration(
      x = csc_score[[j]],
      times = i,
      bandwidth = 0.03, #LOOK AT OTHER METHOD FOR SELECTING BANDWIDTH
      cens.method = "jackknife",
      pseudo = T,
      round = FALSE, 
      xlim = c(0, axislim), ylim = c(0, axislim), 
      rug = TRUE, 
      xlab = "Estimated risks", ylab = "Observed risks",
      auc.in.legend = F, brier.in.legend = F)
    title(paste0("Imputation ",j, ", ", i, "year graft survival (death as competing event)"))
  }
  dev.off()
}



