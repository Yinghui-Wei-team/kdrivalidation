################################################
###              KDRI validation             ###
###       Sample size calculation setup      ###
###              1-year horizon              ###
################################################

#Using the methods of Riley et al. (2022) doi: 10.1002/sim.9275
#Initial code from https://github.com/gscollins1973/Validation_Survival_Sample_Size/blob/main/pesudo_example.R

#Packages
#library(ggplot2)
library(simsurv)
library(survival)


# Time points -------------------------------------------------------------
#This script is run for 1-year graft survival


# # Anticipated LP ----------------------------------------------------------
# #Reported C-index 0.62, which corresponds to Royston's D 0.678
# #Estimate SD using D
# KDRI_sd<- 0.678 / sqrt(8 / pi) # ~0.42487
# #Median KDRI = 1.05
# 
# #Using histogram reported in Rao et al., KDRI can be estimated using LogNormal(mu = log(1.05), 0.42487) 
# set.seed(369)
# KDRI<- rlnorm(10000, meanlog = log(1.05), sdlog = KDRI_sd)
# 
# #Check sampled KDRI looks like histogram reported in Rao et al.
# ggplot(data.frame(KDRI), aes(x=KDRI)) +
#   geom_histogram(aes(y = 100*(..count..)/sum(..count..)),
#                  breaks = seq(0, 5, by = 0.2),
#                  col = "white") + 
#   scale_x_continuous(breaks = seq(0,5,by=0.4)) +
#   labs(y = "Percent",
#        x = "KDRI")


# Survival probabilities --------------------------------------------------
#Using the Kaplan-Meier curves reported in Rao et al. we explore a range of survival probabilities
#Min and max read off plot using DigitizeIt
#1 year
surv_1yr<- seq(0.8745, 0.953, length.out = 4)


# Distribution of event times ---------------------------------------------
#Assume event times follow Exponential distribution
#Find the rate parameter with trial and error
rate_parameter <- function(lambda, LP, time.point){
  X.sim <- simsurv(dist = 'exponential', lambdas = lambda, x = LP, betas = c(KDRI = 1)) #Set beta=1 for calibration slope=1
  X.sim$eventtime <- X.sim$eventtime
  X.sim$dead <- rep(1, length(LP$KDRI))
  X.sim$dead[X.sim$eventtime > time.point] <- 0
  X.sim$eventtime[X.sim$eventtime > time.point] <- time.point
  X.sim$LP <- LP$KDRI

  fit.survfit <- survfit(Surv(eventtime, dead)~1, data = X.sim)
  St          <- summary(fit.survfit, times = time.point)$surv
  fit.coxph   <- coxph(Surv(eventtime, dead)~LP, data = X.sim)
  c.index     <- as.numeric(fit.coxph$concordance[6])
  D           <- as.numeric(royston(fit.coxph)[1])
  return(list(St = St, c.index = c.index, D = D))
}

#1-year graft survival
set.seed(369)
rate_parameter(lambda = 0.010, LP = data.frame(KDRI), time.point = 1)
#Rate 0.045 -> S(1)=0.8558; C=0.64; D=0.90
#Rate 0.010 -> S(1)=0.9634; C=0.66; D=0.95
#This engulfs range of survival proabilities of interest
lambda_1yr<- seq(0.010, 0.045, length.out = 20)

#SEND TO CLUSTER
#Small simulation to account for sampling variability -> gives us lambda
#Initialise
N.SIM<- 100
N<- 5000
St<- matrix(ncol = length(lambda_1yr), nrow = N.SIM)
c.index<- matrix(ncol = length(lambda_1yr), nrow = N.SIM)
D<- matrix(ncol = length(lambda_1yr), nrow = N.SIM)


#1-year survival
for(j in 1:length(lambda_1yr)){
  for(i in 1:N.SIM){
    KDRI<- rlnorm(N, meanlog = log(1.05), sdlog = KDRI_sd)
    OUT<- rate_parameter(lambda = lambda_1yr[j], LP = data.frame(KDRI), time.point = 1)
    St[i,j]<- OUT$St
    c.index[i,j]<- OUT$c.index
    D[i,j]<- OUT$D
  }
}
apply(St, 2, mean)
apply(c.index, 2, mean)
apply(D, 2, mean)


#Get lambdas for range of survival probabilities
St_Interpolate<- approx(x = colMeans(St), y = lambda_1yr, xout = surv_1yr)
C_Interpolate<- approx(x = lambda_1yr, y = colMeans(c.index), xout = approx(x = colMeans(St), y = lambda_1yr, xout = surv_1yr)$y)
D_Interpolate<- approx(x = lambda_1yr, y = colMeans(D), xout = approx(x = colMeans(St), y = lambda_1yr, xout = surv_1yr)$y)

#Collect to write out
samplesurv_1yr<- data.frame(lambda = St_Interpolate$x,
                            surv = St_Interpolate$y,
                            c = C_Interpolate$y,
                            d = D_Interpolate$y)















# #"49,691 patients (72% of total sample) had a functioning graft at end of follow-up"
# #Work with 72% censoring rate
# #Find rate of censoring distribution using method of Wan (2016) doi: 10.1002/sim.7178
# #To avoid estimating the density of the LP, we use the mean LP as it yields the same results
# censor.prop<- function(rate.cens, prop, meanLP) {
#   cen.P<- rate.cens / (rate.cens + meanLP)
#   return(cen.P-prop)
# }
# 
# cens.prop<- 0.72
# mean_LP<- mean(LP)
# 
# rate.cens<- uniroot(censor.prop, prop=cens.prop, meanLP=mean_LP, c(0.001,50), tol=1e-7)
# lambda_cens<- rate.cens$root
# #Rate 2.93800
# 
# 
# 
# 
# 
# 
# #Little simulation to check censoring rate and survival probabilities
# for(j in 1:length(lambda.grid)){  ## This is slow
#   for(i in 1:N.SIM){
#     pb$tick()
#     LP           <- data.frame(id=1:N, value = PearsonDS::rpearson(N, moments = moments))
#     OUT          <- rate_parameter(lambda = lambda.grid[j], LP = LP, time.point = 3)
#     St[i,j]      <- OUT$St
#     c.index[i,j] <- OUT$c.index
#     D[i,j]       <- OUT$D
#   }
# }
# apply(St, 2, mean)
# apply(c.index, 2, mean)
# apply(D, 2, mean)
# 









