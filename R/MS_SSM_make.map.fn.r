make.map.fn = function(data)
{

  #### Prepare options for the map argument of the objective function ####
  map = list()
  map$log_sd_log_NAA<- matrix(1:((data$max_A-1)*data$sp), nrow=data$max_A-1, ncol=data$sp,byrow=T)
  map$log_N1<- array(1:(data$max_A*data$sp), dim = c(data$max_A, data$sp))
  map$log_NAA<- array(1:((data$Y-1)*(data$max_A-1)*data$sp), dim = c(data$Y-1,data$max_A-1,data$sp))
  map$logit_s_surv<-array(1:(data$max_A*data$sp*data$n_surv_max),dim=c(data$max_A,data$sp,data$n_surv_max))
  map$logit_s_F<-array(1:(data$max_A*data$sp),dim=c(data$max_A,data$sp))
  map$logit_q<- array(1:(data$sp*data$n_surv_max), dim = c(data$sp, data$n_surv_max))
  map$logit_gamma_surv<- array(1:(data$sp*data$n_surv_max), dim = c(data$sp, data$n_surv_max))
  map$logit_A50_surv<- array(1:(data$sp*data$n_surv_max), dim = c(data$sp, data$n_surv_max))
  map$acomp_pars_temp_catch<- array(1:(data$sp*3), dim = c(data$sp, 3)) #parameter necessary for age comp distributions, max=3 per species
  map$acomp_pars_temp_index<- array(1:(data$sp*3*data$n_surv_max), dim = c(data$sp,3,data$n_surv_max)) #parameter necessary for age comp distributions, max=3 per species
  map$log_sd_log_rec<-1:data$sp
  map$logit_gamma_F<- 1:data$sp
  map$logit_A50_F<- 1:data$sp
  
  # Reduce the number of parameter to estimate for variance process error survival, CAN BE CHANGED!!!!!!!!!!
  for (i in 2:nrow(map$log_sd_log_NAA))
    map$log_sd_log_NAA[i,] <- map$log_sd_log_NAA[1,] # same variance for all ages


  ## can't estimate random effects and fixed effects for numbers at age in age classes greater than the plus group for that species.
  for(i in 1:data$sp){
    if(data$Aplus[i]<data$max_A){
      map$log_N1[(data$Aplus[i]+1):data$max_A,i] <- NA
      map$log_NAA[,(data$Aplus[i]):(data$max_A-1),i] <- NA
      map$log_sd_log_NAA[(data$Aplus[i]):(data$max_A-1),i] <- NA
    }
    # cannot estimate selectivity after max age catch and survey
    if(data$max_A_catch[i]<data$max_A) map$logit_s_F[(data$max_A_catch[i]+1):data$max_A,i] <- NA
    for (k in 1:data$n_surv[i]){
      if(data$max_A_surv[i,k]<data$max_A) {
        map$logit_s_surv[(data$max_A_surv[i,k]+1):data$max_A,i,k] <- NA
      }
    }
    # cannot estimate selectivity before min age catch and survey
    if(data$min_A_catch[i]>1) map$logit_s_F[1:(data$min_A_catch[i]-1),i] <- NA
    if(data$min_A_catch[i]==data$max_A_catch[i]){ # if only 1 age caught then need to fix selectivity to 1
      #init$logit_s_F[data$min_A_catch[i],i]<-1000
      map$logit_s_F[data$min_A_catch[i],i] <- NA
    }
    for (k in 1:data$n_surv[i]){
      if(data$min_A_surv[i,k]>1) {
        map$logit_s_surv[1:(data$min_A_surv[i,k]-1),i,k] <- NA
      }
      if(data$min_A_surv[i,k]==data$max_A_surv[i,k]){ # if only 1 age caught then need to fix acomp_pars parameters
        #init$logit_s_surv[data$min_A_surv[i,k],i,k]<-1000
        map$acomp_pars_temp_index[i,,k] <- NA
      }
    }
  }
  ## can't estimate survey parameters when survey not available for a species
  for(i in 1:data$sp){
    if(data$n_surv[i]<data$n_surv_max){
      map$logit_q[i,(data$n_surv[i]+1):data$n_surv_max] <- NA 
      map$logit_gamma_surv[i,(data$n_surv[i]+1):data$n_surv_max] <- NA 
      map$logit_A50_surv[i,(data$n_surv[i]+1):data$n_surv_max] <- NA
      map$acomp_pars_temp_index[i,,(data$n_surv[i]+1):data$n_surv_max] <-NA
      map$logit_s_surv[,i,(data$n_surv[i]+1):data$n_surv_max] <- NA
    }
  }
  ## Different number of parameters to estimate for age comp distribution
  for (i in 1:data$sp){
    if (data$age_comp_model_catch[i]==1) # acomp_pars not estimated and should be 0 in init so exp(0)=1
      map$acomp_pars_temp_catch[i,] <- NA # no parameter
    if (data$age_comp_model_catch[i]==2 | data$age_comp_model_catch[i]==3 | data$age_comp_model_catch[i]==5){
      map$acomp_pars_temp_catch[i,2:3] <- NA # only 1 parameter
    }
    if (data$age_comp_model_catch[i]==6){
      map$acomp_pars_temp_catch[i,3] <- NA  # 2 parameters
    }
  }
  for (i in 1:data$sp){
    if (data$age_comp_model_indices[i]==1) # acomp_pars not estimated and should be 0 in init so exp(0)=1
      map$acomp_pars_temp_index[i,,] <- NA # no parameters
    if(data$age_comp_model_indices[i]==2 | data$age_comp_model_indices[i]==3 | data$age_comp_model_indices[i]==5){
      map$acomp_pars_temp_index[i,2:3,] <- NA # only 1 parameter
    }
    if (data$age_comp_model_indices[i]==6){
      map$acomp_pars_temp_index[i,3,] <- NA  # 2 parameters
    }
  }
  ## Different options for M
  map$log_M<-array(1:(data$Y*data$max_A*data$sp),dim=c(data$Y,data$max_A,data$sp))
  map$log_M1<-array(1:(data$max_A*data$sp),dim=c(data$max_A,data$sp))
  map$log_MAA<-array(1:((data$Y-1)*data$max_A*data$sp),dim=c((data$Y-1),data$max_A,data$sp))
  map$log_lorenzen1<-1:data$sp
  map$lorenzen2<-1:data$sp
  map$log_sd_log_MAA<-1:data$sp
  map$logit_scale_M<-array(1:(data$max_A*data$sp),dim=c(data$max_A,data$sp))
  for (i in 1:data$sp){
    #Cannot estimate above the Age + group
    if(data$Aplus[i]<data$max_A){
      map$log_M[,(data$Aplus[i]+1):data$max_A,i] <- NA
      map$log_M1[(data$Aplus[i]+1):data$max_A,i] <- NA
      map$log_MAA[,(data$Aplus[i]+1):data$max_A,i] <- NA
      map$logit_scale_M[(data$Aplus[i]+1):data$max_A,i] <- NA
    }
    #Following depends on M_model
    if (data$M_model[i]==1 || data$M_model[i]==4 || data$M_model[i]==5){
      map$log_M[,,i]<-NA
      map$log_M1[,i]<-NA
      map$log_MAA[,,i]<-NA
      map$log_lorenzen1[i]<-NA
      map$lorenzen2[i]<-NA
      map$log_sd_log_MAA[i]<-NA
    }
    if (data$M_model[i]==2){
      if (data$process_M==1){
        map$log_M[,,i]<-NA
        map$log_lorenzen1[i]<-NA
        map$lorenzen2[i]<-NA
      } else {
        map$log_M1[,i]<-NA
        map$log_MAA[,,i]<-NA
        map$log_lorenzen1[i]<-NA
        map$lorenzen2[i]<-NA
        map$log_sd_log_MAA[i]<-NA
      }
    }
    if (data$M_model[i]==3){
      map$log_M[,,i]<-NA
      map$log_M1[,i]<-NA
      map$log_MAA[,,i]<-NA
      map$log_sd_log_MAA[i]<-NA
    }
    if (data$M_model[i]!=4){
      map$logit_scale_M[,i]<-NA
    }
    # if want the scale on M in M_model=4 to be the same for all ages
    if (data$M_model[i]==4){
      for (i in 1:data$sp)
        map$logit_scale_M[2:data$Aplus[i],i]<-map$logit_scale_M[1,i]
    }
    # Cannot estimate age comp parameters when no age comp data for all years
    for (k in 1:data$n_surv_max){
      if (sum(data$flag_Iaa[,i,k])==0){ 
        map$acomp_pars_temp_index[i,,k]<-NA
      }
    }
  }

  # Dealing with selectivity models
  for (k in 1:data$n_surv_max){
    for (i in 1:data$sp){
      if (data$sel_model_surv[i,k]==1){
        map$logit_A50_surv[i,k]<-NA
        map$logit_gamma_surv[i,k]<-NA
        # #deciding on number of ages estimated, can be changed
        # x_surv=data$max_A_surv[i,k] #default x=data$max_A_surv[i,k] because last selectivity is
        # #x_surv=3
        # for (a in x_surv:data$max_A){
        #   map$logit_s_surv[a,i,k]<-NA
        # }
      }
      if (data$sel_model_surv[i,k]==2){
        map$logit_s_surv[,i,k]<-NA
      }
      if (data$sel_model_F[i]==1){
        map$logit_A50_F[i]<-NA
        map$logit_gamma_F[i]<-NA
        #deciding on number of ages estimated, can be changed
        # x_F=data$Aplus[i] #default x=data$Aplus[i]
        # #x=3
        # for (a in x_F:data$max_A){
        #   map$logit_s_F[a,i]<-NA
        # }
        
      }
      if (data$sel_model_F[i]==2){
        map$logit_s_F[,i]<-NA
      }
    }
  }
  
  ## Deal with the process errors
  if (data$process_survival==0){
    map$log_NAA[]<-NA
    map$log_sd_log_NAA[]<-NA
  }
  if (data$process_rec==0){
    map$log_sd_log_rec[]<-NA
  }
  if(data$process_rec==1 & data$recruit_model==2){
    map$mean_log_rec<-1:data$sp
  } else {
    map$mean_log_rec<-rep(NA,data$sp)
  }
  for (i in 1:data$sp){
    if (data$M_model[i]==2 & data$process_M==0){
      # for (j in 2:data$Y){
      #   map$log_M[j,,i]<-map$log_M[1,,i] # M cst over time but estimated across ages
      # }
      for (a in 2:data$Aplus[i]){
        map$log_M[,a,i]<-map$log_M[,1,i] # M cst over age but estimated over time
      }
    }
    if (data$M_model[i]==2 & data$process_M==1){
      for (i in 1:data$sp)
          map$log_M1[,i]<-1
    }
  }
  
  
  # For trophic interactions
  if (data$predation_on==0){ # if not predation all interaction parameters = NAs
    map$vuln_par <- matrix(nrow=data$sp,ncol=data$n_pred)
    map$log_scale_gamma_par <- rep(NA,data$n_pred)
    map$log_shape_gamma_par <- rep(NA,data$n_pred)
    map$log_power_typeIII <- rep(NA,data$n_pred)
    map$par_deltadir <- array(NA,dim=c(data$Y,3,data$n_pred))
    map$log_cons_rate <- array(NA,dim=c(data$max_B,data$n_pred))
    map$log_sd_cons_rate <- array(NA,dim=c(data$max_B,data$n_pred))
    map$log_alpha_cons <- array(NA, dim=c(data$sp,data$n_pred))
    map$log_beta_cons <- array(NA, dim=c(data$sp,data$n_pred))
  } else { # if predation is on, depends on assumptions
    map$vuln_par <-matrix(1:(data$sp*data$n_pred),nrow=data$sp,ncol=data$n_pred)
    for (i in 1:data$sp){
      for (j in 1:data$n_pred){
        if (data$interaction[i,j]==0){
          map$vuln_par[i,j]<-NA #vuln not estimated if no interaction between predator and prey
        }
      }
    }
    if (data$gamma_pref_estim==0){
      map$log_scale_gamma_par<-rep(NA,data$n_pred)
      map$log_shape_gamma_par<-rep(NA,data$n_pred)
    }
    if (data$functional_response==1){
      map$log_power_typeIII<-rep(NA,data$n_pred) # no estimated power when type II functional response
    } else {
      if (data$functional_response!=2){
        stop("functional_response option does not exist")
      }
    }
    if (data$cons_rate_estim==0){
      map$log_cons_rate <- array(NA,dim=c(data$max_B,data$n_pred))
      map$log_sd_cons_rate <- array(NA,dim=c(data$max_B,data$n_pred))
      map$log_alpha_cons <- array(NA, dim=c(data$sp,data$n_pred))
      map$log_beta_cons <- array(NA, dim=c(data$sp,data$n_pred))
    } else {
      map$log_cons_rate <- array(1:(data$max_B*data$n_pred),dim=c(data$max_B,data$n_pred))
      map$log_sd_cons_rate <- array(1:(data$max_B*data$n_pred),dim=c(data$max_B,data$n_pred))
      # #if want to assume same consumption rate for all prey
      # for(i in 2:data$sp){
      #   map$log_cons_rate[i,,] <- map$log_cons_rate[1,,]
      #   map$log_sd_cons_rate[i,,] <- map$log_sd_cons_rate[1,,]
      # }
      #if want to assume same consumption rate for all ages
      for(b in 2:data$max_B){
        map$log_cons_rate[b,] <- map$log_cons_rate[1,]
        map$log_sd_cons_rate[b,] <- map$log_sd_cons_rate[1,]
      }
      ## can't estimate for age classes greater than the plus group for that predator
      for(j in 1:data$n_pred) if(data$Bplus[j]<data$max_B){ 
        map$log_cons_rate[(data$Bplus[j]+1):data$max_B,j] <- NA
        map$log_sd_cons_rate[(data$Bplus[j]+1):data$max_B,j] <- NA
      }
      map$log_alpha_cons <- array(1:(data$sp*data$n_pred), dim=c(data$sp,data$n_pred))
      map$log_beta_cons <- array(1:(data$sp*data$n_pred), dim=c(data$sp,data$n_pred))
      ## if want to assume same alpha and beta for all prey
      if (data$sp>1) {
        for(i in 2:data$sp){
        map$log_alpha_cons[i,] <- map$log_alpha_cons[1,]
        map$log_beta_cons[i,] <- map$log_beta_cons[1,]
        # map$log_cons_rate[i,,] <- map$log_cons_rate[1,,]
        # map$log_sd_cons_rate[i,,] <- map$log_sd_cons_rate[1,,]
        }
      }
    }
    
    map$par_deltadir <- array(1:(data$Y*3*data$n_pred),dim=c(data$Y,3,data$n_pred))
    # if want to assume same delta-dirichlet for all years for phi and intercept (i in 1:3), just phi (i in 1:1)
    for (t in 2:data$Y){
      map$par_deltadir[t,,]=map$par_deltadir[1,,]
    }
    # if want to fix slope parameter of delta-dirichlet
    map$par_deltadir[,2,] <- NA
    #map$par_deltadir[,,] <- NA
    
    #if want to fix the general vulnerability
    #map$vuln_par <- matrix(nrow=data$sp,ncol=data$n_pred)
    
  }
  if (data$biomass_other_option==1){
    map$log_biomass_other_y1<-NA
    map$log_K_other<-NA
    map$log_growth_other<-NA
  }
  
  # ### Not estimating cons_rate_other for now
  # map$log_cons_rate_other<-array(dim=c(data$max_B,data$n_pred))
  # map$log_sd_cons_rate_other<-array(dim=c(data$max_B,data$n_pred))
  
  if (data$diet_model==2) map$par_deltadir[,2:3,]=NA

  if (data$process_F==0){
    map$logit_phi_process_F <- rep(NA,data$sp)
    map$log_sd_process_F <- rep(NA,data$sp)
    map$cst_process_F <- rep(NA,data$sp)
  }
  
  ## Cannot estimate q if selectivity estimated at age (because confounded) and q=1
  for (i in 1:data$sp){
    for (k in 1:data$n_surv[i]){
      if (data$sel_model_surv[i,k]==1){
        map$logit_q[i,k]=NA
      }
    }
  }
  
  ## Need to fix one s_F if random F and s_F estimated by age
  if (data$process_F==1){
    for (i in 1:data$sp) {
      if(data$sel_model_F[i]==1){ # need to fix one value so selectivity and random effect not confounded
        map$logit_s_F[data$max_A_catch[i],i]<- NA # selectivity of 1 in last catch age
      }
    }
  }
  
  ## option if burn-in TRUE
  if (data$burnin==1) map$log_N1[]=NA # N1 not estimated

  
  #### test
  map$logit_s_F[data$max_A_catch-1,]=NA
  map$logit_s_F[data$max_A_catch,]=NA
  # map$logit_scale_M[]=1:20
  # map$log_M[,6:10,1]=map$log_M[,5,1]
  # map$par_deltadir[]=NA
  # map$log_sd_log_NAA[]=NA
  # map$vuln_par=matrix(NA,nrow=data$sp,ncol=data$n_pred)
  
  #map$logit_q[]=NA
  #map$logit_A50_surv[]=NA
  #map$logit_gamma_surv[]=NA
  #map$acomp_pars_temp_index[]=NA
  #map$logit_s_surv[]=NA
  #map$log_N1[20]=NA
  #map$acomp_pars_temp_catch[]=NA

  #map$par_deltadir[,3,]=NA
  #map$logit_s_surv[1,1,]=NA
  #map$mean_log_rec[]=NA
  #map$lorenzen2[]=NA
  #map$log_lorenzen1[]=NA
  #map$logit_s_surv[2,2,2]<-NA
  # map$log_alpha_cons[]=1
  # map$log_beta_cons[]=1
  # map$log_cons_rate[]=1
  # map$log_sd_cons_rate[]=1
  #map$log_biomass_other_y1<-NA
  #map$log_K_other<-NA
  #map$log_growth_other<-NA
  # map$log_E<-array(1:(data$Y*data$sp),dim=c(data$Y,data$sp))
  # map$log_rec<-array(1:(data$Y-1*data$sp),dim=c(data$Y-1,data$sp))
  # map$log_E[40:42,2]<-NA
  # map$log_sd_log_NAA[]<-NA
  # map$log_sd_log_rec[]<-NA
  ####
  map = lapply(map, function(x) factor(x))
  return(map)
  
}

