make.random.fn = function(data)
{

  if (data$process_rec==1 & data$process_survival==1){
    random<-c("log_rec","log_NAA")
  }
  if (data$process_rec==1 & data$process_survival==0){
    random<-c("log_rec")
  }
  if (data$process_rec==0 & data$process_survival==1){
    random<-c("log_NAA")
  }
  if (data$process_rec==0 & data$process_survival==0){
    random<-NULL
  }
  if (data$process_M==1){
    if (data$M_model[1]==2)
      random<-c(random,"log_MAA")
    if (data$M_model[1]==3)
      random<-c(random,"log_lorenzen1","lorenzen2")
  }
  
  # if (data$gamma_pref_estim == 1){
  #   random<-c(random,"log_shape_gamma_par","log_scale_gamma_par")
  # }
  # 
  if (data$predation_on==1){
    if (data$cons_rate_estim == 1){
    random<-c(random,"log_alpha_cons","log_beta_cons")
    }
  }
  
  if (data$process_F==1){
    random<-c(random,"log_E")
  }
  
  return(random)
}

