#### Estimate Convertion parameters to convert length to age ####

load("data/NOAA/LENGTH_W_DATA_FROM_SOLE_10SP.RDATA")

all.sp <- c(73,72,32,74,197,15,23,105,121,106) # SVSPP in order "cod", "silverhake", "herring", "haddock", "goosefish", "spinydogfish","winterskate", "yellowtailflounder", "mackerel","winterflounder"
sp=all.sp[1:data$n_pred] # only extracted for number of predators


data1=bio.data[which(bio.data$SVSPP==sp),]
data2<-data1[which(is.na(data1$AGE)==FALSE),] # remove NAs in age
data2<-data2[which(data2$AGE!=0),] # remove 0s in age so can take logs


lm_coeff <- matrix(nrow=2,ncol=data$n_pred)
for (j in 1:data$n_pred){
  data_pred <- data2[which(data2$SVSPP==sp[j]),]
  reg <- lm(log(data_pred$AGE)~log(data_pred$LENGTH))
  lm_coeff[,j] <- reg$coefficients
}
# reg <- lm(log(data2$AGE)~log(data2$LENGTH))
# plot(log(data_pred$AGE)~log(data_pred$LENGTH))
# abline(reg)
# print(summary(reg))
# 
# plot(data_pred$AGE~data_pred$LENGTH)
# length=1:200
# age=exp(log(length)*reg$coefficients[2]+reg$coefficients[1])
# lines(age~length,col="red")
# lines(log(age)~log(length),col="red")


# #### Determine age for each stomach ####
# 
data$age_pred <- diet_data$length_pred

for (t in 1:data$Y){
  for (j in 1:data$n_pred) {
    for (n in 1:data$n_stom[t,j]) {
      data$age_pred[n,t,j]<-round(exp(log(diet_data$length_pred[n,t,j])*lm_coeff[2,j]+lm_coeff[1,j]))
      if (data$age_pred[n,t,j]==0) data$age_pred[n,t,j] <- 1 # change age 0 into age 1
    }
  }
}

