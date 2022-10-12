###Sample size

library(pmsampsize)

#Data
data<- read.csv("../../NHSBT Data/Data/kdri_mi.csv", na.strings = "")
D<- filter(data, X_mi_m %in% 0:10)
original<- filter(data, X_mi_m==0)

#Time horizon
timehor<- c(1,5)

#Number of candidate predictor parameters
p<- 13

#Mean follow-up
meanfup<- mean(original$tsurv/365.5, na.rm = T)

#Event rate
E<- c("1year" = filter(original, gsurv<365.5 & gcens==1) %>% nrow(), #1236
      "5years" = filter(original, gsurv<(365.5*5) & gcens==1) %>% nrow()) #2392

eventrate<- c("1year" = 1236/20035, "5year" = 2392/20035)

#R^2 using C-stat reported by Watson
#This process is detailed in Riley et al. (2020) Min sample size paper
C<- 0.63
D<- 5.5*(C-0.5) + (10.26*(C-0.5)^3)
R2D<- (pi/8*D^2) / ((pi^2/6) + (pi/8*D^2))
R2oQ<- ((-1) * pi^2/6 * R2D) / (((1 - (pi^2)/6) * R2D) - 1)
LR<- (-1) * E * log(1 - R2oQ)
R2CSapp<- 1 - exp(-LR/E)
SVH<- 1 - (p/LR) #p candidate predictor parameters
R2CSadj<- SVH * R2CSapp 



#Sample size calculation
sampsize<- data.frame()
for(i in 1:length(timehor)) {
  samp<- pmsampsize(type = "s",
                    rsquared = R2CSadj[i],
                    parameters = p,
                    rate = eventrate[i],
                    timepoint = timehor[i],
                    meanfup = meanfup)
  
  sampsize<- rbind(sampsize, data.frame(time = timehor[i],
                                        eventrate = eventrate[i],
                                        R2 = samp$r2a,
                                        sample_size = samp$sample_size,
                                        events = ceiling(samp$events),
                                        EPP = samp$EPP))
}

write.csv(sampsize, "../output/samplesize_calc.csv", row.names = F)

