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

V <- rep(0,3529000)
t1 <- 1
dim(V) <- c(1000,3529)
N <- rep(0,3529000)
t2 <- 1
dim(N) <- c(1000,3529)
for (i in 1:2000){
  summfitID <- as.numeric(summary(survfit(myCPH, newdata = data.frame(x1= x1[i], x2= x2[i], x3=x3[i])))$surv)
  if (d[i]==1){
     V[t1,1:3529]<-summfitID[1:3529]
     t1 <- t1+1
  }else{
     N[t2,1:3529]<-summfitID[1:3529]
     t2 <- t2+1
  }
}

hist(V[1:t1,100],100)
hist(N[1:t2,200],100)
summfitID <- as.numeric(summary(survfit(myCPH, newdata = data.frame(x1= 10, x2= 0.9, x3=2)))$surv)


