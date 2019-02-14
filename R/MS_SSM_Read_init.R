read.ssm_init.fn <- function(initf)
{
  
  char.lines <- readLines(initf) # read inita file line by line
  com.ind <- which(substring(char.lines,1,1) == "#") # line number of comment's lines
  init.start <- com.ind[c(which(diff(com.ind)>1), length(com.ind))] # Find comment line where inita follows and remove others
  #comments <- char.lines[init.start] # entire line of comments
  
  # Create a list of all the observation inita by extracting them from the inita file
  init <- list()
  ind <- 0
  
  gamma_F <- scan(initf, skip = init.start[ind <- ind + 1], n = data$sp)
  init$logit_gamma_F <- -log((data$max_A_catch/gamma_F)-1)
  
  A50_F <- scan(initf,  skip = init.start[ind <- ind + 1], n = data$sp)
  init$logit_A50_F <- -log((data$max_A_catch/A50_F)-1)
  
  N1<-matrix(scan(initf, what = double(), skip = init.start[ind <- ind + 1], n = data$max_A*data$sp), data$max_A, data$sp, byrow = TRUE)
  init$log_N1 <- log(N1)
  
  rec<-matrix(scan(initf, what = double(), skip = init.start[ind <- ind + 1], n = (data$Y-1)*data$sp), (data$Y-1), data$sp, byrow = TRUE)
  init$log_rec <- log(rec)
  
  E<-matrix(scan(initf, what = double(), skip = init.start[ind <- ind + 1], n = data$Y*data$sp), data$Y, data$sp, byrow = TRUE)
  init$log_E <- log(E)
  
  NAA.temp<-lapply(1:data$sp, function(x) matrix(scan(initf, what = double(), skip = init.start[ind+x], n = (data$Y-1)*(data$max_A-1)), data$Y-1, data$max_A-1, byrow = TRUE))
  init$log_NAA<-array(dim=c(data$Y-1,data$max_A-1,data$sp))
  for (i in 1:data$sp)
    init$log_NAA[,,i]<-log(NAA.temp[[i]])
  ind <- ind+data$sp
  
  cv_log_NAA<-matrix(scan(initf, what = double(), skip = init.start[ind <- ind + 1], n = (data$max_A-1)*data$sp), data$max_A-1, data$sp, byrow = TRUE)
  init$log_sd_log_NAA<-log(sqrt(log(cv_log_NAA^2+1)))
  
  cv_log_rec<-scan(initf, skip = init.start[ind <- ind + 1], n = data$sp)
  init$log_sd_log_rec<-log(sqrt(log(cv_log_rec^2+1)))
  
  q <- matrix(scan(initf, what = double(), skip = init.start[ind <- ind + 1], n = data$sp*data$n_surv_max), data$sp, data$n_surv_max, byrow = TRUE)
  ## Cannot estimate q if selectivity estimated at age (because confounded)
  for (i in 1:data$sp){
    for (k in 1:data$n_surv[i]){
      if (data$sel_model_surv[i,k]==1){
        q[i,k]=0.999999999
      }
    }
  }
  init$logit_q <- log(q/(1-q))
  
  
  gamma_surv <- matrix(scan(initf, what = double(), skip = init.start[ind <- ind + 1], n = data$sp*data$n_surv_max), data$sp, data$n_surv_max, byrow = TRUE)
  init$logit_gamma_surv <- -log((data$max_A/gamma_surv)-1)
  
  A50_surv <- matrix(scan(initf, what = double(), skip = init.start[ind <- ind + 1], n = data$sp*data$n_surv_max), data$sp, data$n_surv_max, byrow = TRUE)
  init$logit_A50_surv <- -log((data$max_A/A50_surv)-1)
  
  agecomp_pars<- matrix(scan(initf, what = double(), skip = init.start[ind <- ind + 1], n = 6*3), 6, 3, byrow = TRUE) #6=total number of possible distributions for age composition, 3=max number of parameters for distribution
  init$acomp_pars_temp_catch <- rbind(agecomp_pars[data$age_comp_model_catch,]) #(3)
  acomp_pars_temp_index <- rbind(agecomp_pars[data$age_comp_model_indices,]) #(3)
  init$acomp_pars_temp_index<-array(dim=c(data$sp,3,data$n_surv_max))
  for (i in 1:data$n_surv_max)
    init$acomp_pars_temp_index[,,i]<-acomp_pars_temp_index
  
  mean_rec <- scan(initf, skip = init.start[ind <- ind + 1], n = data$sp)
  init$mean_log_rec <- log(mean_rec)
  
  #### For M options, no need to modify the init file ##########
  init$log_M<-log(data$M)
  
  init$log_M1<-cbind(log(data$M[1,,]))
  
  init$log_MAA<-array(dim=c(data$Y-1,data$max_A,data$sp))
  for (i in 1:data$sp)
    init$log_MAA[,,i]<-log(data$M[-1,,i])
  
  init$log_lorenzen1<-rep(log(3.69),data$sp)
  
  init$lorenzen2<-rep(-0.305,data$sp)
  ##############################################################
  
  cv_log_MAA<-scan(initf, skip = init.start[ind <- ind + 1], n = data$sp)
  init$log_sd_log_MAA<-log(sqrt(log(cv_log_MAA^2+1)))
  
  scale_M <- matrix(scan(initf, what = double(), skip = init.start[ind <- ind + 1], n = data$sp*data$max_A), data$max_A, data$sp, byrow = TRUE)
  init$logit_scale_M <- log(scale_M/(data$scale_M_upper-scale_M))
  
  init$logit_s_surv <- array(0,dim=c(data$max_A,data$sp,data$n_surv_max)) # selectivity=1 by default
  # for (i in 1:data$sp){
  #   for (k in 1:data$n_surv[i]){
  #     if(data$min_A_surv[i,k]==data$max_A_surv[i,k]){ # if only 1 age caught then need to fix selectivity to 1
  #       init$logit_s_surv[data$min_A_surv[i,k],i,k]<-1000
  #       #map$logit_s_surv[data$min_A_surv[i,k],i,k] <- NA
  #     }
  #   }
  # }
  
  init$logit_s_F <- array(0,dim=c(data$max_A,data$sp)) # selectivity=1 by default
  for (i in 1:data$sp){
    if(data$min_A_catch[i]==data$max_A_catch[i]){ # if only 1 age caught then need to fix selectivity to 1
      init$logit_s_F[data$min_A_catch[i],i]<-1000
      #map$logit_s_F[data$min_A_catch[i],i] <- NA
    }
  }
  
  
  vuln_par<-rep(1/data$n_prey,data$sp) 
  init$vuln_par <- array(dim=c(data$sp,data$n_pred))
  for (j in 1:data$n_pred)
    init$vuln_par[,j] <- vuln_par # vulnerability=for all prey and predators by default

  init$log_shape_gamma_par <- data$log_shape_gamma_dat
  
  init$log_scale_gamma_par <- data$log_scale_gamma_dat
  
  init$log_power_typeIII <- rep(0,data$n_pred) # by default power_typeIII=1 so type II
  
  deltadir<-scan(initf, skip = init.start[ind <- ind + 1], n = 3)
  init$par_deltadir<-array(dim=c(data$Y,3,data$n_pred))
  for (i in 1:data$Y)
    init$par_deltadir[i,,]<-deltadir
  
  init$log_cons_rate<-matrix(scan(initf, what = double(), skip = init.start[ind <- ind + 1], n = data$max_B*data$n_pred), data$max_B, data$n_pred, byrow = TRUE)
  # init$log_cons_rate <- array(dim=c(data$max_B,data$n_pred))
  # for (j in 1:data$n_pred)
  #   for (i in 1:data$sp)
  #     init$log_cons_rate[i,,j]<-log_cons_rate[i,j]
  
  #log_sd_cons_rate<-matrix(scan(initf, what = double(), skip = init.start[ind <- ind + 1], n = data$sp*data$n_pred), data$sp, data$n_pred, byrow = TRUE)
  # init$log_sd_cons_rate <- array(dim=c(data$sp,data$max_B,data$n_pred))
  # for (j in 1:data$n_pred)
  #   for (i in 1:data$sp)
  #     init$log_sd_cons_rate[i,,j]<-log_sd_cons_rate[i,j]
  log_sd_cons_rate<-scan(initf, what = double(), skip = init.start[ind <- ind + 1], n =data$n_pred)
  init$log_sd_cons_rate <- array(dim=c(data$max_B,data$n_pred))
  for (j in 1:data$n_pred)
    init$log_sd_cons_rate[,j]<-log_sd_cons_rate[j]
  
  
  alpha_cons <- scan(initf, skip = init.start[ind <- ind + 1], n = 1)
  beta_cons <- scan(initf, skip = init.start[ind <- ind + 1], n = 1)
  init$log_alpha_cons<-matrix(log(alpha_cons),nrow=data$sp,ncol=data$n_pred)
  init$log_beta_cons<-matrix(log(beta_cons),nrow=data$sp,ncol=data$n_pred)
  
  
  ### other food
  # log_cons_rate_other<-scan(initf, skip = init.start[ind <- ind + 1], n = data$n_pred)
  # init$log_cons_rate_other <- array(dim=c(data$max_B,data$n_pred))
  # for (i in 1:data$n_pred)
  #   init$log_cons_rate_other[,i]<-log_cons_rate_other[i]
  # 
  # log_sd_cons_rate_other<-scan(initf, skip = init.start[ind <- ind + 1], n = data$n_pred)
  # init$log_sd_cons_rate_other <- array(dim=c(data$max_B,data$n_pred))
  # for (i in 1:data$n_pred)
  #   init$log_sd_cons_rate_other[,i]<-log_sd_cons_rate_other[i]
  
  biomass_other_y1<-scan(initf, skip = init.start[ind <- ind + 1], n = 1)
  init$log_biomass_other_y1 <- log(biomass_other_y1)
  
  growth_other<-scan(initf, skip = init.start[ind <- ind + 1], n = 1)
  init$log_growth_other <- log(growth_other)
  
  K_other<-scan(initf, skip = init.start[ind <- ind + 1], n = 1)
  init$log_K_other <- log(K_other)
  ###
  
  sd_process_F<-scan(initf, skip = init.start[ind <- ind + 1], n = data$sp)
  init$log_sd_process_F <- log(sd_process_F)
  
  phi_process_F<-scan(initf, skip = init.start[ind <- ind + 1], n = data$sp)
  init$logit_phi_process_F <- -log(2/(phi_process_F+1)-1) # comprised between -1 and 1
  
  mean_E<-scan(initf, skip = init.start[ind <- ind + 1], n = data$sp)
  init$cst_process_F <- log(mean_E)*(1-phi_process_F)
  
  
  return(init)
  
}