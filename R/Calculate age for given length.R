#### Estimate Convertion parameters to convert length to age ####

load("data/NOAA/LENGTH_W_DATA_FROM_SOLE_10SP.RDATA")

sp=73 # cod


data1=bio.data[which(bio.data$SVSPP==sp),]
data2<-data1[which(is.na(data1$AGE)==FALSE),] # remove NAs in age
data2<-data2[which(data2$AGE!=0),] # remove 0s in age so can take logs


reg <- lm(log(data2$AGE)~log(data2$LENGTH))
# plot(log(data2$AGE)~log(data2$LENGTH))
# abline(reg)
# print(summary(reg))
# 
# plot(data2$AGE~data2$LENGTH)
# length=1:200
# age=exp(log(length)*reg$coefficients[2]+reg$coefficients[1])
# lines(age~length,col="red")
# lines(log(age)~log(length),col="red")


lm_coeff=reg$coefficients


# #### Determine age for each stomach ####
# 
data$age_pred <- diet_data$length_pred

for (t in 1:data$Y){
  for (n in 1:data$n_stom[t,1]) { # 1 because just for cod
    data$age_pred[n,t,1]<-round(exp(log(diet_data$length_pred[n,t,1])*lm_coeff[2]+lm_coeff[1]))
    if (data$age_pred[n,t,1]==0) data$age_pred[n,t,1] <- 1 # change age 0 into age 1
  }
}

