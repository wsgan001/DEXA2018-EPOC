prediction <- rep(-1, nrow(index_small))
for ( i in 1:nrow(index_small)){
  index <- index_small[i,]
  data = time_small[as.numeric(index[2]):as.numeric(index[3]),]
  friend = friend_small[as.numeric(index[2]):as.numeric(index[3]),]
  temp <- length(data)/2
  pred.time <- seq(0, 600000, by = 60)
  infectiousness <- get.infectiousness(data[1:temp], friend[1:temp], pred.time)
  pred <- pred.cascade(pred.time, infectiousness$infectiousness, data[1:temp], friend[1:temp])
  #prediction[i] <- pred[length(pred)]*2
  prediction[i] <- abs(index_small[i,3]-index_small[i,2]-2*pred[length(pred)]-1)*500+pred[length(pred)]
  print (i)
}


data(tweet)
pred.time <- seq(0, 6 * 60 * 60, by = 60)
infectiousness <- get.infectiousness(tweet[, 1], tweet[, 2], pred.time)
plot(pred.time, infectiousness$infectiousness)
data(tweet)
pred.time <- seq(0, 6 * 60 * 60, by = 60)
infectiousness <- get.infectiousness(tweet[, 1], tweet[, 2], pred.time)
pred <- pred.cascade(pred.time, infectiousness$infectiousness, tweet[, 1], tweet[, 2], n.star = 100)
plot(pred.time, pred)
