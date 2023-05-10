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


# Time points -------------------------------------------------------------
#This script is run for 5-year graft survival


# Anticipated LP ----------------------------------------------------------
#Reported C-index 0.62, which corresponds to Royston's D 0.678
#Estimate SD using D
KDRI_sd<- 0.678 / sqrt(8 / pi) # ~0.42487
#Median KDRI = 1.05

#Using histogram reported in Rao et al., KDRI can be estimated using LogNormal(mu = log(1.05), 0.42487) 
set.seed(369)
KDRI<- rlnorm(10000, meanlog = log(1.05), sdlog = KDRI_sd)

#Check sampled KDRI looks like histogram reported in Rao et al.
ggplot(data.frame(KDRI), aes(x=KDRI)) +
  geom_histogram(aes(y = 100*(..count..)/sum(..count..)),
                 breaks = seq(0, 5, by = 0.2),
                 col = "white") + 
  scale_x_continuous(breaks = seq(0,5,by=0.4)) +
  labs(y = "Percent",
       x = "KDRI")


# Survival probabilities --------------------------------------------------
#Using the Kaplan-Meier curves reported in Rao et al. we explore a range of survival probabilities
#Min and max read off plot using DigitizeIt
surv_5yr<- seq(0.635, 0.822, length.out = 4)



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

#5-year graft survival
rate_parameter(lambda = 0.010, LP = data.frame(KDRI), time.point = 5)
#Rate 0.028 -> S(5)=0.6304; C=0.64; D=0.90
#Rate 0.010 -> S(5)=0.8294; C=0.63; D=0.80
lambda_5yr<- seq(0.010,0.028, length.out = 20)

#Small simulation to account for sampling variability -> gives us lambda
#Initialise
N.SIM<- 10
N<- 5000
St<- matrix(ncol = length(lambda_1yr), nrow = N.SIM)
c.index<- matrix(ncol = length(lambda_1yr), nrow = N.SIM)
D<- matrix(ncol = length(lambda_1yr), nrow = N.SIM)

pb <- progress::progress_bar$new(format = "  simulation :what [:bar] :percent eta: :eta",
                                 clear = FALSE, 
                                 total = N.SIM * length(lambda_1yr), 
                                 width = 60)

#1-year survival
for(j in 1:length(lambda_1yr)){
  for(i in 1:N.SIM){
    pb$tick()
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

plot(lambda_1yr, colMeans(St), type = 'b', pch = 20, xlab = 'lambda_c', ylab = 'S(3)')
grid()

#Get lambdas for range of survival probabilities
approx(x = colMeans(St), y = lambda_1yr, xout = surv_1yr)
approx(x = lambda_1yr, y = colMeans(c.index), xout = approx(x = colMeans(St), y = lambda_1yr, xout = surv_1yr)$y)
approx(x = lambda_1yr, y = colMeans(D), xout = approx(x = colMeans(St), y = lambda_1yr, xout = surv_1yr)$y)














