# extended cox 
#install.packages("survival")
#install.packages("OIsurv")
library(OIsurv)
N <- dim(relapse)[1]
t1 <- rep(0, N+sum(!is.na(relapse$int))) # initialize start time at 0
t2 <- rep(-1, length(t1)) # build vector for end times
d <- rep(-1, length(t1)) # whether event was censored
g <- rep(-1, length(t1)) # gender covariate
i <- rep(FALSE, length(t1)) # initialize intervention at FALSE

j <- 1
for(ii in 1:dim(relapse)[1]){
  if(is.na(relapse$int[ii])){ # no intervention, copy survival record
    t2[j] <- relapse$event[ii]
    d[j] <- relapse$delta[ii]
    g[j] <- relapse$gender[ii]
    j <- j+1
  } else { # intervention, split records
    g[j+0:1] <- relapse$gender[ii] # gender is same for each time
    d[j] <- 0 # no relapse observed pre-intervention
    d[j+1] <- relapse$delta[ii] # relapse occur post-intervention?
    i[j+1] <- TRUE # intervention covariate, post-intervention
    t2[j] <- relapse$int[ii]-1 # end of pre-intervention
    t1[j+1] <- relapse$int[ii]-1 # start of post-intervention
    t2[j+1] <- relapse$event[ii] # end of post-intervention
    j <- j+2 # two records added
  }
}

mySurv <- Surv(t1, t2, d) # pg 3 discusses left-trunc. right-cens. data
myCPH <- coxph(mySurv ~ g + i)
summfitID <- as.numeric(summary(survfit(myCPH, newdata = data.frame(g = 0,i = TRUE)))$surv)
summfitID
plot(summfitID)

t1 <- ok$V4
t2 <- ok$V5
d <- ok$V6
x1 <- ok$V1
x2 <- ok$V2
x3 <- ok$V3
mySurv <- Surv(t1, t2, d) # pg 3 discusses left-trunc. right-cens. data
myCPH <- coxph(mySurv ~ x1+x2+x3)
summary(myCPH)
summfitID <- as.numeric(summary(survfit(myCPH, newdata = data.frame(x1= 10, x2= 0.9, x3=2)))$surv)
summfitID <- as.numeric(summary(survfit(myCPH))$surv)
summfitID
plot(summfitID)


result <- rep(-1, 200)
for(i in 1:200){
  summfitID <- as.numeric(summary(survfit(myCPH, newdata = data.frame(x1= ok[i,1], x2= ok[i,2], x3= ok[i,3])))$surv)
  result[i] <- summfitID[20]
}
hist(result,20)
