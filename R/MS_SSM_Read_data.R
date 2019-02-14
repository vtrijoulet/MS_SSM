read.ssm_dat.fn <- function(datf)
{
  
  char.lines <- readLines(datf) # read data file line by line
  com.ind <- which(substring(char.lines,1,1) == "#") # line number of comment's lines
  dat.start <- com.ind[c(which(diff(com.ind)>1), length(com.ind))] # Find comment line where data follows and remove others
  #comments <- char.lines[dat.start] # entire line of comments
  
  # Create a list of all the observation data by extracting them from the data file
  dat <- list()
  ind <- 0
  
  
  dat$sp <- scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = 1) #number of modelled species
  dat$n_prey<-dat$sp+1 # modelled species + other food
  
  #dat$sp_names <- scan(datf, what = character(), skip = dat.start[ind <- ind + 1], n = dat$sp)
  
  dat$interaction <- matrix(scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$sp*dat$sp), dat$sp, dat$sp, byrow = TRUE)
  
  dat$n_pred <- length(which(apply(dat$interaction,2,sum)>0)) #number of predator species

  dat$year1 <- scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$sp) #first year with catch data per species
  
  dat$lastyear <- scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = 1) #last year for data (same for all sp, possibility to change it to a vector for different last year per species)

  dat$Y <- dat$lastyear-min(dat$year1)+1 #position of last year (possibility to change it to a vector for different last year per species)
  dat$Y1<- dat$year1-min(dat$year1)+1 #position first year per species

  dat$Aplus <- scan(datf, what = integer(), skip = dat.start[ind <- ind + 1],n=dat$sp) #age in age+ group for preys
  
  dat$Bplus <- dat$Aplus[which(apply(dat$interaction,2,sum)>0)] #age in age+ group for predators
  
  dat$max_A<-max(dat$Aplus) #max number of age classes for prey for array dimensions
  dat$max_B<-max(dat$Bplus) #max number of age classes for predators for array dimensions
  
  dat$min_A_catch <-  scan(datf, what = integer(), skip = dat.start[ind <- ind + 1],n=dat$sp) 
  dat$max_A_catch <-  scan(datf, what = integer(), skip = dat.start[ind <- ind + 1],n=dat$sp) 
  
  dat$min_Fbar_idx <-  scan(datf, what = integer(), skip = dat.start[ind <- ind + 1],n=dat$sp) 
  dat$max_Fbar_idx <-  scan(datf, what = integer(), skip = dat.start[ind <- ind + 1],n=dat$sp) 
  dat$Fbar_range <- dat$max_Fbar_idx-dat$min_Fbar_idx+1
  
  month_spawn <- scan(datf, what = double(), skip = dat.start[ind <- ind + 1], n =dat$sp)
  dat$prop_y_elapsed_SSB<-(month_spawn-1)/12
  
  wc.temp<-lapply(1:dat$sp, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$max_A), dat$Y, dat$max_A, byrow = TRUE))
  dat$wc<-array(dim=c(dat$Y,dat$max_A,dat$sp))
  for (i in 1:dat$sp)
    dat$wc[,,i]<-wc.temp[[i]]
  ind <- ind+dat$sp
  
  wssb.temp<-lapply(1:dat$sp, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$max_A), dat$Y, dat$max_A, byrow = TRUE))
  dat$wSSB<-array(dim=c(dat$Y,dat$max_A,dat$sp))
  for (i in 1:dat$sp)
    dat$wSSB[,,i]<-wssb.temp[[i]]
  ind <- ind+dat$sp
  
  mat.temp<-lapply(1:dat$sp, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$max_A), dat$Y, dat$max_A, byrow = TRUE))
  dat$mature<-array(dim=c(dat$Y,dat$max_A,dat$sp))
  for (i in 1:dat$sp)
    dat$mature[,,i]<-mat.temp[[i]]
  ind <- ind+dat$sp
  
  dat$obs_aggr_Cw<-matrix(scan(datf, what = double(), skip = dat.start[ind <- ind + 1], n = dat$Y*dat$sp), dat$Y, dat$sp, byrow = TRUE)
  cv_obs_aggr_Cw<-matrix(scan(datf, what = double(), skip = dat.start[ind <- ind + 1], n = dat$Y*dat$sp), dat$Y, dat$sp, byrow = TRUE)
  #cv_obs_aggr_Cw<-cv_obs_aggr_Cw+0.5 #test increase in CV for catch
  dat$sd_obs_aggr_Cw<-sqrt(log(cv_obs_aggr_Cw^2+1))
  
  dat$flag_aggr_Cw<-matrix(scan(datf, what = double(), skip = dat.start[ind <- ind + 1], n = dat$Y*dat$sp), dat$Y, dat$sp, byrow = TRUE)
  
  M.temp<-lapply(1:dat$sp, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$max_A), dat$Y, dat$max_A, byrow = TRUE))
  dat$M<-array(dim=c(dat$Y,dat$max_A,dat$sp))
  for (i in 1:dat$sp)
       dat$M[,,i]<-M.temp[[i]]
  ind <- ind+dat$sp
  
  dat$n_surv_max <- scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = 1) #number of different surveys
  
  dat$n_surv <- scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$sp) #number of surveys per species
  
  dat$min_A_surv <- matrix(scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$sp*dat$n_surv_max), dat$sp, dat$n_surv_max, byrow = TRUE)
  dat$max_A_surv <- matrix(scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$sp*dat$n_surv_max), dat$sp, dat$n_surv_max, byrow = TRUE)
  
  month_surv<-matrix(scan(datf, what = double(), skip = dat.start[ind <- ind + 1], n = dat$Y*dat$n_surv_max), dat$Y, dat$n_surv_max, byrow = TRUE)
  dat$prop_y_elapsed_surv<-(month_surv-1)/12
  
  I.temp<-lapply(1:dat$n_surv_max, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$sp), dat$Y, dat$sp, byrow = TRUE))
  dat$obs_aggr_I<-array(dim=c(dat$Y,dat$sp,dat$n_surv_max))
  for (i in 1:dat$n_surv_max)
    dat$obs_aggr_I[,,i]<-I.temp[[i]]
  ind <- ind+dat$n_surv_max

  CV.temp<-lapply(1:dat$n_surv_max, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$sp), dat$Y, dat$sp, byrow = TRUE))
  CV.temp2<-array(dim=c(dat$Y,dat$sp,dat$n_surv_max))
  for (i in 1:dat$n_surv_max)
    CV.temp2[,,i]<-CV.temp[[i]]
  #CV.temp2[,,2]<-CV.temp2[,,2]+0.5 #test increase in CV for surveys
  dat$sd_obs_aggr_I<-sqrt(log(CV.temp2^2+1))
  ind <- ind+dat$n_surv_max
  
  flag.temp<-lapply(1:dat$n_surv_max, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$sp), dat$Y, dat$sp, byrow = TRUE))
  dat$flag_aggr_I<-array(dim=c(dat$Y,dat$sp,dat$n_surv_max))
  for (i in 1:dat$n_surv_max)
    dat$flag_aggr_I[,,i]<-flag.temp[[i]]
  ind <- ind+dat$n_surv_max
  
  propCaa.temp<-lapply(1:dat$sp, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$max_A), dat$Y, dat$max_A, byrow = TRUE))
  dat$obs_prop_Caa<-array(dim=c(dat$Y,dat$max_A,dat$sp))
  for (i in 1:dat$sp)
    dat$obs_prop_Caa[,,i]<-propCaa.temp[[i]]
  ind <- ind+dat$sp
  
  dat$flag_Caa<-matrix(scan(datf, what = double(), skip = dat.start[ind <- ind + 1], n = dat$Y*dat$sp), dat$Y, dat$sp, byrow = TRUE)  

  propIaa_surv1.temp<-lapply(1:dat$sp, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$max_A), dat$Y, dat$max_A, byrow = TRUE))
  if(dat$n_surv_max>1){
    ind <- ind+dat$sp
    propIaa_surv2.temp<-lapply(1:dat$sp, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$max_A), dat$Y, dat$max_A, byrow = TRUE))
  }
  if(dat$n_surv_max>2){
    ind <- ind+dat$sp
    propIaa_surv3.temp<-lapply(1:dat$sp, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$max_A), dat$Y, dat$max_A, byrow = TRUE))
  }
  if(dat$n_surv_max>3){
    ind <- ind+dat$sp
    propIaa_surv4.temp<-lapply(1:dat$sp, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$max_A), dat$Y, dat$max_A, byrow = TRUE))
  }
  if(dat$n_surv_max>4){
    ind <- ind+dat$sp
    propIaa_surv5.temp<-lapply(1:dat$sp, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$max_A), dat$Y, dat$max_A, byrow = TRUE))
  }
  if(dat$n_surv_max>5){
    ind <- ind+dat$sp
    propIaa_surv6.temp<-lapply(1:dat$sp, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$max_A), dat$Y, dat$max_A, byrow = TRUE))
  }
  
  dat$obs_prop_Iaa<-array(dim=c(dat$Y,dat$max_A,dat$sp,dat$n_surv_max))
  for (i in 1:dat$sp){
    dat$obs_prop_Iaa[,,i,1]<-propIaa_surv1.temp[[i]]
    if(dat$n_surv_max>1){
      dat$obs_prop_Iaa[,,i,2]<-propIaa_surv2.temp[[i]]
    }
    if(dat$n_surv_max>2){
      dat$obs_prop_Iaa[,,i,3]<-propIaa_surv3.temp[[i]]
    }
    if(dat$n_surv_max>3){
      dat$obs_prop_Iaa[,,i,4]<-propIaa_surv4.temp[[i]]
    }
    if(dat$n_surv_max>4){
      dat$obs_prop_Iaa[,,i,5]<-propIaa_surv5.temp[[i]]
    }
    if(dat$n_surv_max>5){
      dat$obs_prop_Iaa[,,i,6]<-propIaa_surv6.temp[[i]]
    }
  }
  ind <- ind+dat$sp
  
  flag.temp2<-lapply(1:dat$n_surv_max, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$sp), dat$Y, dat$sp, byrow = TRUE))
  dat$flag_Iaa<-array(dim=c(dat$Y,dat$sp,dat$n_surv_max))
  for (i in 1:dat$n_surv_max)
    dat$flag_Iaa[,,i]<-flag.temp2[[i]]
  ind <- ind+dat$n_surv_max
  
  dat$Neff_C<-matrix(scan(datf, what = double(), skip = dat.start[ind <- ind + 1], n = dat$Y*dat$sp), dat$Y, dat$sp, byrow = TRUE)
  
  #dat$Neff_surv<-matrix(scan(datf, what = double(), skip = dat.start[ind <- ind + 1], n = dat$Y*dat$n_surv_max), dat$Y, dat$n_surv_max, byrow = TRUE)
  Neff.temp<-lapply(1:dat$n_surv_max, function(x) matrix(scan(datf, what = double(), skip = dat.start[ind+x], n = dat$Y*dat$sp), dat$Y, dat$sp, byrow = TRUE))
  dat$Neff_surv<-array(dim=c(dat$Y,dat$sp,dat$n_surv_max))
  for (i in 1:dat$n_surv_max)
    dat$Neff_surv[,,i]<-Neff.temp[[i]]
  ind <- ind+dat$n_surv_max
  
  dat$age_comp_model_catch <- scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$sp) #index for choice of LL distribution for age comp catch
  
  dat$age_comp_model_indices <- scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$sp) #index for choice of LL distribution for age comp survey indices
  
  dat$scale_M_upper<-scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$sp)
  
  dat$sel_model_surv<-matrix(scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$sp*dat$n_surv_max), dat$sp, dat$n_surv_max, byrow = TRUE)
  
  dat$sel_model_F<-scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = dat$sp)
  
  dat$log_shape_gamma_dat<-scan(datf, skip = dat.start[ind <- ind + 1], n = dat$n_pred)
  
  dat$log_scale_gamma_dat<-scan(datf, skip = dat.start[ind <- ind + 1], n = dat$n_pred)
  
  dat$log_ratio_w_diet<-as.matrix(read.table("data/NOAA/log_ratio_w_diet_cod.txt"))
  dat$n_ratio<-nrow(dat$log_ratio_w_diet)
  
  dat$B_ecosystem<-7340000 #tons
  
  dat$B_other<-scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = 1)
  
  ny_aggr<-scan(datf, what = integer(), skip = dat.start[ind <- ind + 1], n = 1) # number of years for aggregation diet flag_nll_diet==4
  n_aggr <- dat$Y/ny_aggr
  diff <- n_aggr-floor(n_aggr)
  # If last interval >=5 years keep it otherwise increase size last interval:
  if (diff<0.5) {
    n_aggr<-round(n_aggr)
    dat$t1 <- seq(0,dat$Y,ny_aggr)[1:n_aggr]
    dat$t2 <- c(seq(ny_aggr,dat$Y,ny_aggr)[1:(n_aggr-1)],dat$Y)
    } else {
    n_aggr<-ceiling(n_aggr)
    dat$t1 <- seq(0,dat$Y,ny_aggr)
    dat$t2 <- c(seq(ny_aggr,dat$Y,ny_aggr),dat$Y)
  }
  dat$length_t<-dat$t2-dat$t1 
  dat$xmax<-length(dat$t1)


  return(dat)

}
