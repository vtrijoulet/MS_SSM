data$catch_aref = matrix(data$max_A, data$Y, data$sp)
data$index_aref = array(data$max_A, dim=c(data$Y, data$sp,data$n_surv_max))
# data$Neff_C[]=500
# data$Neff_surv[]=500
# data$age_comp_model_catch[]=1
# data$age_comp_model_indices[]=1
# for (i in 1:data$sp){
#   if (data$age_comp_model_catch[i]==1){
#     init$acomp_pars_temp_catch[i,1]=0 #Needs to be 0 so exp(0)=1
#   }
#   if (data$age_comp_model_indices[i]==1){
#     init$acomp_pars_temp_index[i,1,]=0
#   }
# }
#init$par_deltadir[,3,]=20
data$ratio_diet3=array(1,dim=c(data$Y,data$sp+1,data$max_B,data$n_pred))
#init$par_deltadir[,1,]=5 #################### 5=variance almost 0, perfect diet data
vuln=rbind(0.450,0.309,0.24)
vuln=cbind(vuln,rbind(NA,0.50,0.49)) # original
#vuln=cbind(vuln,rbind(NA,0.01,0.989)) # new vuln
#vuln=rbind(0.759,0.24) # when only cod as predator
vuln.other=c(0.001,0.01) # original
#vuln.other=c(0.001,0.001) # new vuln
sum.vuln=1/vuln.other-1
init$vuln_par[]=log(vuln*(1+sum.vuln))
source('R/Calculate_age_for_given_length.R')
#data$sd_obs_aggr_I[]<-sqrt(log(0.3^2+1))



# data$obs_aggr_Cw[] = exp(6)
# data$obs_aggr_I[] = exp(6)
data$spring_cons<-array(dim=c(data$Y_cons,2,data$n_pred))
data$fall_cons<-array(dim=c(data$Y_cons,2,data$n_pred))
data$spring_cons[,,1]<-spring[,,1,1]
data$fall_cons[,,1]<-fall[,,1,1]
data$spring_cons_other<-array(dim=c(data$Y_cons,2,data$n_pred))
data$fall_cons_other<-array(dim=c(data$Y_cons,2,data$n_pred))
data$spring_cons_other[,,1]<-spring.other[,,1]
data$fall_cons_other[,,1]<-fall.other[,,1]
data$flag_cons_rate<-rep(1,data$Y_cons)
for(t in 1:data$Y_cons){
  for (i in 1:data$sp){
    if (sum(is.na(data$spring_cons[t,,1]))>0 || sum(is.na(data$fall_cons[t,,1]))>0)
      data$flag_cons_rate[t]<-0 
  }
}



# upper.par=vector(length=length(obj2$par))
# upper.par[]=Inf
# upper.par[which(names(obj2$par)=="mean_log_rec")]=20
# system.time(opt2<-nlminb(obj2$par, obj2$fn, obj2$gr, upper=upper.par, control = list(eval.max = 5000, iter.max = 5000))) #,sing.tol=1e-20)))
# 
# 
# 
# 
# 
# init$log_N1[,1]=seq(12,7.5,by=-0.5)
# init$log_N1[,2]=seq(16,11.5,by=-0.5)
# init$log_rec[,1]=12
# init$log_rec[,2]=16
# init$log_NAA[,,1]=10
# init$log_NAA[,,2]=14
# init$mean_log_rec[1]=12
# init$mean_log_rec[1]=16
# data$B_other=150000
# diet_data<-extract_diet_data_fn("~/NOAA/Foodhabitsdata/Brian's data/diet_GB_All_10sp_with_empty.csv")
# data$ratio_diet <- array(dim=c(575,data$Y,(data$n_prey),data$n_pred)) #modelled species sp + other food
# data$prob_l_given_b <- diet_data$prob_l_given_b # (l,Y,max_B,n_pred)
# data$n_stom=array(575,dim=c(data$Y,data$n_pred))
# data$n_stom_max <- max(data$n_stom)
# data$flag_cons_rate<-array(1,dim=c(data$Y_cons))
# data$spring_cons<-array(dim=c(data$Y_cons,2,data$sp,data$n_pred))
# data$fall_cons<-array(dim=c(data$Y_cons,2,data$sp,data$n_pred))
# data$spring_cons_other<-array(dim=c(data$Y_cons,2,data$n_pred))
# data$fall_cons_other<-array(dim=c(data$Y_cons,2,data$n_pred))
# 
# 
# 
# source("MS_SSM_extract_diet_data.R") 
# diet_data<-extract_diet_data_fn("~/NOAA/Foodhabitsdata/Brian's data/diet_GB_All_10sp_with_empty.csv")
# data$prob_l_given_b <- diet_data$prob_l_given_b # (l,Y,max_B,n_pred)
# source('~/NOAA/Multispecies SSM/MS_SSM_test_NLL_in_R.R')
# 
# ratio_diet<-array(0,dim=dim(data$prob_l_given_b))
