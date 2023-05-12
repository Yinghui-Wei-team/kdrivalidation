################################################
###              KDRI validation             ###
###       Sample size calculation setup      ###
###              5-year horizon              ###
################################################

#Using the methods of Riley et al. (2022) doi: 10.1002/sim.9275
#Initial code from https://github.com/gscollins1973/Validation_Survival_Sample_Size/blob/main/pesudo_example.R

#Packages
library(ggplot2)
library(simsurv)
library(survival)



# Save location -----------------------------------------------------------
#Computer
# folder<- "samplesize/"
#Cluster
folder<- "kdrival/"


# Time points -------------------------------------------------------------
#This script is run for 5-year graft survival


# Anticipated LP ----------------------------------------------------------
#Reported C-index 0.62, which corresponds to Royston's D 0.678
#Estimate SD using D
KDRI_sd<- 0.678 / sqrt(8 / pi) # ~0.42487
#Median KDRI = 1.05

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
# 
# LP<- data.frame(id = 1:N.OBS,
#                 value = log(KDRI))


# Survival probabilities --------------------------------------------------
#Using the Kaplan-Meier curves reported in Rao et al. we explore a range of survival probabilities
#Min and max read off plot using DigitizeIt
surv_5yr<- seq(0.635, 0.822, length.out = 4)



# Distribution of event times ---------------------------------------------
#Assume event times follow Exponential distribution
#Find the rate parameter with trial and error
rate_parameter <- function(lambda, LP, time.point){
  X.sim <- simsurv(dist = 'exponential', lambdas = lambda, x = LP, betas = c(value = 1)) #Set beta=1 for calibration slope=1
  X.sim$eventtime <- X.sim$eventtime
  X.sim$dead <- rep(1, length(LP$value))
  X.sim$dead[X.sim$eventtime > time.point] <- 0
  X.sim$eventtime[X.sim$eventtime > time.point] <- time.point
  X.sim$LP <- LP$value
  
  fit.survfit <- survfit(Surv(eventtime, dead)~1, data = X.sim)
  St          <- summary(fit.survfit, times = time.point)$surv
  fit.coxph   <- coxph(Surv(eventtime, dead)~LP, data = X.sim)
  c.index     <- as.numeric(fit.coxph$concordance[6])
  D           <- as.numeric(royston(fit.coxph)[1])
  return(list(St = St, c.index = c.index, D = D))
}


#5-year graft survival
# set.seed(369)
# rate_parameter(lambda = 0.032, LP = LP, time.point = 5)
#Rate 0.089 -> S(5)~0.6105; C=0.62; D=0.69
#Rate 0.032 -> S(5)~0.8369; C=0.62; D=0.72
lambda_5yr<- seq(0.032,0.089, length.out = 20)

#Small simulation to account for sampling variability -> gives us lambda
#Initialise
set.seed(369)
N.SIM<- 100
N.OBS<- 10000
St<- matrix(ncol = length(lambda_5yr), nrow = N.SIM)
c.index<- matrix(ncol = length(lambda_5yr), nrow = N.SIM)
D<- matrix(ncol = length(lambda_5yr), nrow = N.SIM)

# pb <- progress::progress_bar$new(format = "  simulation [:bar] :percent eta: :eta",
#                                  clear = FALSE, 
#                                  total = N.SIM * length(lambda_5yr), 
#                                  width = 60)

#5-year survival
for(j in 1:length(lambda_5yr)){
  for(i in 1:N.SIM){
    #pb$tick()
    if (i %in% seq(10,N.SIM, length.out=10)) {print(paste0("sim ", i))}
    
    KDRI<- rlnorm(N.OBS, meanlog = log(1.05), sdlog = KDRI_sd)
    LP<- data.frame(id = 1:N.OBS,
                    value = log(KDRI))
    OUT<- rate_parameter(lambda = lambda_5yr[j], LP = LP, time.point = 5)
    St[i,j]<- OUT$St
    c.index[i,j]<- OUT$c.index
    D[i,j]<- OUT$D
  }
}
apply(St, 2, mean)
apply(c.index, 2, mean)
apply(D, 2, mean)

# plot(lambda_5yr, colMeans(St), type = "b", pch = 20, xlab = "lambda", ylab = "S(5)")
# grid()

#Get lambdas for range of survival probabilities
St_Interpolate<- approx(x = colMeans(St), y = lambda_5yr, xout = surv_5yr)
C_Interpolate<- approx(x = lambda_5yr, y = colMeans(c.index), xout = approx(x = colMeans(St), y = lambda_5yr, xout = surv_5yr)$y)
D_Interpolate<- approx(x = lambda_5yr, y = colMeans(D), xout = approx(x = colMeans(St), y = lambda_5yr, xout = surv_5yr)$y)

#Collect to write out
samplesurv_5yr<- data.frame(surv = St_Interpolate$x,
                            lambda = St_Interpolate$y,
                            c = C_Interpolate$y,
                            d = D_Interpolate$y)
write.csv(samplesurv_5yr, paste0(folder, "lambdas_5yr.csv"), row.names = F)


# Distribution of censoring times -----------------------------------------
#"49,691 patients (72% of total sample) had a functioning graft at end of follow-up"
#Work with 72% censoring rate
#Find rate of censoring distribution using method of Wan (2016) doi: 10.1002/sim.7178
#To avoid estimating the density of the LP, we use the mean LP as it yields the same results
censor.prop<- function(rate.cens, prop, meanLP) {
  cen.P<- rate.cens / (rate.cens + meanLP)
  return(cen.P-prop)
}

#Simulation to estimate lambda
print("cens sim")
N.SIM<- 500
N.OBS<- 5000
cens.prop<- 0.72
lambda_cens<- matrix(ncol=1, nrow=N.SIM)
cens<- matrix(ncol=1, nrow=N.SIM)
# pb <- progress_bar$new(format = "  simulation [:bar] :percent eta: :eta",
#                        clear = FALSE, 
#                        total = N.SIM, 
#                        width = 60)


for (i in 1:N.SIM) {
  #pb$tick()
  if (i %in% seq(10,N.SIM, length.out=10)) {print(paste0("cens sim ", i))}
  
  KDRI<- rlnorm(N.OBS, meanlog = log(1.05), sdlog = KDRI_sd)
  mean_LP<- mean(log(KDRI))
  
  #Find the lambda which provides 72% censoring on average
  rate.cens<- uniroot(censor.prop, prop=cens.prop, meanLP=mean_LP, c(0.001,50), tol=1e-7)
  lambda_cens[i,1]<- rate.cens$root
  
  #Double check this gives correct censoring rate
  X.sim<- simsurv(dist = 'exponential', lambdas = lambda_cens[i,1], x = data.frame(ids = seq(1:N.OBS)))
  X.sim$cens<- rep(1, nrow(X.sim))
  X.sim$cens[X.sim$eventtime > 5]<- 0
  X.sim$eventtime<- X.sim$eventtime
  X.sim$eventtime[X.sim$eventtime > 5]<- 5
  fit.survfit<- survfit(Surv(eventtime, status)~1, data = X.sim)
  cens[i,1]<- summary(fit.survfit, times = 4.99999999)$sur
}

mean_censprob<- apply(cens, 2, mean)
mean_lambdacens<- apply(lambda_cens, 2, mean)

print("mean cens prob")
mean_censprob

print("mean lambda for cens")
mean_lambdacens










