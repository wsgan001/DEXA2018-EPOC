# extended cox 
#install.packages("survival")
#install.packages("OIsurv")
library(OIsurv)

t1 <- ok$V4
t2 <- ok$V5
d <- ok$V6
x1 <- ok$V1
x2 <- ok$V2
x3 <- ok$V3
mySurv <- Surv(t1, t2, d) # pg 3 discusses left-trunc. right-cens. data
myCPH <- coxph(mySurv ~ x1+x2+x3)
summary(myCPH)
#summfitID <- as.numeric(summary(survfit(myCPH, newdata = data.frame(x1= 10, x2= 0.9, x3=2)))$surv)
summfitID <- as.numeric(summary(survfit(myCPH))$surv)
summfitID
plot(summfitID)
m <- summary(survfit(myCPH))$cumhaz
write.csv(m,'cumhaz.csv')

