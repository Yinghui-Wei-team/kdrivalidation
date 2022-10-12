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
data<- read.csv("../../NHSBT Data/Data/kdri_mi.csv", na.strings = "")

D<- data %>%
  filter(X_mi_m > 0 & X_mi_m <=10) %>%
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

#dage* dheight* dweight* dethnic* dpast_hypertension dpast_diabetes dcreat* dhcv dtype

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


#Recalibrate the KDRI using the cause-specific Cox model
#Number of imputations
n.imp<- 10

#Model formula
outcome<- "Hist(etime, etype_fct)"
vars<- c("rcdage18", "rcdage", "rcdage50", "rcdheight", "rcdweight",
         "dethnicb_i", "dpast_hypertension_i", "dpast_diabetes_i", "cva_i",
         "rcdcreat", "rcdcreat15", "dhcv_i", "dtype_i")

kdriform<- as.formula(paste(outcome,
                            paste(vars, collapse = "+"),
                            sep = "~"))

#Hold all model output in list
rckdri_csc<- list()

for(i in 1:n.imp) {
  #Data for imputation i
  tempD<- filter(D_cr, X_mi_m==i)
  
  #Cause-specific Cox model with variables used in the KDRI
  rckdri_csc[[i]]<- CSC(kdriform, data = tempD, cause = 1)
}

coefs<- data.frame()
for(i in 1:n.imp) {
  coefs<- rbind(coefs, data.frame(imputation=i,
                                  var = vars,
                                  summary(rckdri_csc[[i]]$models$`Cause 1`)$coefficients[,c("coef", "se(coef)")]))
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

write.csv(pool_coef, "../output/recalibration/by_imputation/csc_poolcoef.csv", row.names = F)

#Linear predictor
D_cr$LP<- as.matrix(D_cr[,vars]) %*% as.matrix(pool_coef[,"pool_coef"])

#Distribution of LP
ggplot(D_cr, aes(x=LP-uskdri_lp)) +
  geom_histogram() +
  facet_wrap(~X_mi_m)


#Validate using a competing risks model
#Hold all model output in list
rckdrival_csc<- list()

#C index
csc_cstat<- data.frame()

B<- 500 #Bootstrap samples
csc_cboot<- data.frame(matrix(nrow=B)) #Data frame to store C-index from each bootstrap sample

#Store results from score function
csc_score<- list()

#Store AUC and Brier scores
csc_auc<- data.frame()
csc_brier<- data.frame()

#Store Aalen-Johansen and predicted risks
csc_risks<- data.frame()



#START LOOP
for (i in 1:n.imp) {
  #Data for imputation i
  tempD<- filter(D_cr, X_mi_m==i)
  
  #Include KDRI LP in a cause-specific Cox model
  rckdrival_csc[[i]]<- CSC(Hist(etime, etype_fct) ~ LP, data = tempD,
                           cause = 1)
  
  print("models ran")
  
  #C-index using IPCW because of Wolbers 2014 (C-index defined in Wolbers 2009)
  #Apparent C index
  csc_tempout<- pec::cindex(rckdrival_csc[[i]],
                            formula = Hist(etime, etype_fct) ~ 1,
                            data = tempD,
                            splitMethod = "BootCv",
                            B = B,
                            keep.matrix = T,
                            verbose=F)
  print("Store cindex")
  
  csc_cstat<- rbind(csc_cstat, data.frame(imputation=i,
                                          appCstat=csc_tempout$AppCindex$CauseSpecificCox,
                                          bootCstat=csc_tempout$BootCvCindex$CauseSpecificCox))
  
  print("cindex stored")
  
  csc_cboot<- cbind(csc_cboot, csc_tempout$BootstrapCrossValCindexMat[[1]])
  
  
  #Time-dependent AUC and Brier score at 1-year and 5-years following transplantation
  csc_score[[i]] <- riskRegression::Score(
    list("kdri_val" = rckdrival_csc[[i]]),
    formula = Hist(etime, etype_fct) ~ 1, 
    cause = 1,
    cens.model = "km", 
    data = tempD, 
    conf.int = TRUE, 
    times = c(1,5),
    metrics = c("auc", "brier"),
    summary = c("risks"),
    plots = "calibration",
    split.method = "bootcv",
    B = B)
  
  csc_auc<- rbind(csc_auc, data.frame(imputation=i,
                                      filter(csc_score[[i]]$AUC$score, model=="kdri_val")[,c("times", "AUC", "se")]))
  
  csc_brier<- rbind(csc_brier, data.frame(imputation=i,
                                          filter(csc_score[[i]]$Brier$score, model=="kdri_val")[,c("times", "Brier", "se")]))
  
  
  #Predicted risk from model and pseudo values
  csc_risks<- rbind(csc_risks, data.frame(imputation=i, csc_score[[i]]$Calibration$plotframe))
  
  print(paste0("imputation ",i, " done"))
  
  
}

saveRDS(rckdrival_csc, file = "../output/rc_csc_model.RData")
saveRDS(csc_score, file = "../output/rc_csc_score.RData")


write.csv(csc_auc, file="../output/recalibration/by_imputation/rc_csc_auc.csv", row.names = F)
write.csv(csc_brier, file="../output/recalibration/by_imputation/rc_csc_brier.csv", row.names = F)
write.csv(csc_cstat, file="../output/recalibration/by_imputation/rc_csc_cstat.csv", row.names = F)
write.csv(csc_cboot, file="../output/recalibration/by_imputation/rc_csc_cboot.csv", row.names = F)
write.csv(csc_risks, file="../output/rc_csc_risks.csv", row.names = F)


#Calibration
n.imp<- 10
csc_cal<- data.frame(matrix(nrow = n.imp, ncol = 2))
csc_calboot<- data.frame()
B<- 500

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
    cal_slope<- geepack::geese(pseudovalue ~ offset(clog_risks) + clog_risks,
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
      bootdata<- pseudos_sorted[sample(nrow(pseudos_sorted), replace = TRUE),]
      tempcal_slope<- geepack::geese(pseudovalue ~ offset(clog_risks) + clog_risks,
                                     data = bootdata, id = ID,
                                     scale.fix = TRUE,
                                     family = gaussian,
                                     mean.link = "cloglog",
                                     corstr = "independence",
                                     jack = TRUE)
      csc_calboot<- rbind(csc_calboot, data.frame(imputation=i,
                                                  times=timehor[j],
                                                  calslope=summary(tempcal_slope)$mean["clog_risks","estimate"],
                                                  cal_se= summary(tempcal_slope)$mean["clog_risks","san.se"]))
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
