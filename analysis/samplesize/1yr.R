################################################
###              KDRI validation             ###
###       Sample size calculation setup      ###
###              1-year horizon              ###
################################################

#Using the methods of Riley et al. (2022) doi: 10.1002/sim.9275
#Initial code from https://github.com/gscollins1973/Validation_Survival_Sample_Size/blob/main/pesudo_example.R


# Folder setup ------------------------------------------------------------
#Computer
folder<- "samplesize/"
#Cluster
#folder<- "kdrival/"


# Scenarios ---------------------------------------------------------------
#   | Scenario | Time horizon | Survival probability | Event time rate |
#   |----------|--------------|----------------------|-----------------|
#   | 1        | 1            | 0.8745               | 0.118491        |
#   | 2        | 1            | 0.900667             | 0.092057        |
#   | 3        | 1            | 0.926833             | 0.066684        |
#   | 4        | 1            | 0.953                | 0.042119        |
#   | 5        | 5            | 0.635                | 0.082733        |
#   | 6        | 5            | 0.697333             | 0.064967        |
#   | 7        | 5            | 0.759667             | 0.049146        |
#   | 8        | 5            | 0.822                | 0.034688        |

scenario<- 1

#Select one lambda at a time so that they can run in parallel
survlambdas<- read.csv(paste0(folder, "lambdas_1yr.csv"))
lambda<- survlambdas$lambda[1]

#Rate for censoring distribution
#Found in setup step
lambda_cens<- 0.1257687

#Try for variety of sample sizes 
sampsize<- seq(10000, 24000, by=2000)



# Initialise simulation ---------------------------------------------------
set.seed(369)
N.OBS<- 100000
N.SIM<- 100
time.point<- 1
se.slope<- matrix(ncol = length(sampsize), nrow = N.SIM)
S1<- matrix(ncol = length(sampsize), nrow = N.SIM)
KDRI_sd<- 0.678 / sqrt(8 / pi) # ~0.42487
KDRI<- rlnorm(N.OBS, meanlog = log(1.05), sdlog = KDRI_sd)
LP_all<- data.frame(id = 1:N.OBS,
                    value = log(KDRI))


# Simulation --------------------------------------------------------------
for (j in 1:length(sampsize)) {
  for (i in 1:N.SIM) {
    N.VAL<- sampsize[j]
    X.sim<- simsurv(dist = 'exponential', lambdas = 0.029, x = LP_all, betas = c(value = 1))
    X.sim$LP<- LP_all$value
    
    #Censor at max follow-up (10 years)
    X.sim$status <- rep(1, nrow(X.sim))
    X.sim$status[X.sim$eventtime > 10]    <- 0
    X.sim$eventtime[X.sim$eventtime > 10] <- 10
    
    #Generate censoring times
    X.cens <- simsurv(dist = 'exponential', lambdas = lambda_cens, x = data.frame(ids=seq(1:nrow(X.sim))))
    X.sim$status[X.cens$eventtime < X.sim$eventtime] <- 0
    X.sim$eventtime[X.cens$eventtime < X.sim$eventtime] <- X.cens$eventtime[X.cens$eventtime < X.sim$eventtime]
    
    #Set up calibration plots
    #Want to save as tiff
    tiff(paste0(folder, "calplot_scenario_", scenario, ".tif"))
    plot(c(0,0.4), c(0,0.4), type= 'n', xlab = 'predicted risk', ylab = 'observed risk')
    
    for(i in 1:N.SIM){
      X.sim$val <- rep(0, nrow(X.sim))
      #Define validation cohort randomly
      index <- sample(1:nrow(X.sim), size = N.VAL, replace = F)
      X.sim$val[index] <- 1
      #Fit model in those not invalidation cohort
      fit <- coxph(Surv(eventtime, status)~LP, data = X.sim[X.sim$val == 0, ])
      
      S1[i,j]      <- summary(survfit(fit), times = time.point)$surv
      #Predicted risk for those in validation cohort
      pred.val   <- 1 - S1[i,j]^predict(fit, newdata = X.sim[X.sim$val == 1, ], type = 'risk')
      f.val      <- prodlim(Hist(eventtime, status)~1, data = X.sim[X.sim$val == 1, ])
      pseudo.val <- jackknife(f.val, times = time.point)
      
      X <- matrix(ncol=2, nrow = length(pred.val))
      X[,1] <- pred.val
      X[,2] <- pseudo.val
      X <- X[order(X[,1]),]
      #Plot calibration curves for each simulation
      lines(X[, 1], 1 - predict(loess(X[,2]~X[,1])), col = 'grey')
      
      pv.val <- 1 - pseudo.val ## pseudo values are generated for survival (so 1 minus for the event)
      mv.val <- pred.val
      xx.val <- log(-log(1 - mv.val))
      
      # Fit calibration model
      calfit.val.2 <- geepack::geese(pv.val~xx.val, 
                                     jack      = T, 
                                     scale.fix = T, 
                                     id        = 1:length(xx.val), 
                                     family    = gaussian, 
                                     mean.link = 'cloglog', 
                                     corstr    = 'independence')
      
      se.slope[i,j] <- summary(calfit.val.2)$mean[2,2] # pull out robust se
    }
    abline(a = 0, b = 1, col = 'black', lwd = 2)
    dev.off()
    
  }
}

summary(se.slope)

pdf(file = paste0(folder, "runningcalplot_scenario_", scenario, ".pdf"))
for (j in 1:length(sampsize)) {
  plot(cumsum(se.slope[,j])/seq_along(se.slope[,j]), pch = 20, 
       ylab='running mean', xlab='iteration',
       main = paste0("Sample size ", sampsize[j], "; lambda ", lambda, "; time horizon ", time.point))
}
dev.off()



