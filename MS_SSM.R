######################################################################################################################################
##                                                                                                                                  ##
##                                  Multispecies state-space assessment model                                                       ##
##                                                                                                                                  ##
##      Trijoulet V., Fay G. and Miller T.J. (2019). Performance of a state-space multispecies model:                               ## 
##      what are the consequences of ignoring predation and process errors in stock assessments?                                    ##
##      Journal of Applied Ecology. https://doi.org/10.1111/1365-2664.13515                                                         ##
##                                                                                                                                  ##
##      Trijoulet, V., Fay, G., Curti, K., Smith, B., & Miller, T. J. (2019). Performance of multispecies population models:        ##
##      insights on the influence of diet data. ICES Journal of Marine Science. https://doi.org/10.1093/icesjms/fsz053              ##
##                                                                                                                                  ##
######################################################################################################################################



library ("TMB")

#gdbsource("MS_SSM.R",TRUE) # debugger



#### Compile C++ code and load compiled file ####

#compile("MS_SSM.cpp","-O1 -g",DLLFLAGS="") # if need to debug
compile("src/MS_SSM.cpp")
dyn.load(dynlib("src/MS_SSM"))



#### Source necessary scripts ####

source("R/MS_SSM_Read_data.R")
source("R/MS_SSM_Read_init.R")
source("R/MS_SSM_get_aref_fn.R")
# For trophic interactions
source("R/MS_SSM_extract_diet_data.R")
source("R/MS_SSM_obs_cons_rate.R")
#source("MS_SSM_extract_prob_l_given_b.R")
source("R/MS_SSM_fake_length_age_ratio_diet.R")
source("R/MS_SSM_burnin.fn.R")

#### Read data file ####

data<-read.ssm_dat.fn("data/data_WBSS_herring.dat") # (same names than in C++ code)
init<-read.ssm_init.fn("data/Initial_param_WBSS_herring.dat")

# data<-read.ssm_dat.fn("data/data_cod+silverhake+herring_like.dat") # (same names than in C++ code)
# init<-read.ssm_init.fn("data/Initial_param_cod_silverhake_herring_like.dat")

# Determine aref with get_aref_fn
data$catch_aref = matrix(NA, data$Y, data$sp)
for(i in 1:data$sp)
  data$catch_aref[,i] = get_aref_fn(data$obs_prop_Caa[,1:data$max_A_catch[i],i])

data$index_aref = array(NA, dim=c(data$Y, data$sp,data$n_surv_max))
for(i in 1:data$n_surv_max)
  for(j in 1:data$sp)
    data$index_aref[,j,i] = get_aref_fn(data$obs_prop_Iaa[,1:data$max_A_surv[j,i],j,i])



#### Choose options for model run ####

data$burnin <- 0 # if 0 no burn-in period and log_N1 estimated, contrary if 1

data$recruit_model<-2 # only used when process_rec=1: 1=ramdom walk, 2=random about the mean

data$process_rec<-1 # 1=process error on recruitemnt log_rec, 0= no process error
data$process_survival<-1 # 1=process error on survival log_NAA, 0 = noprocess error
data$process_F <- 1 # 1=AR1 on log_E

data$M_model<-rep(1,data$sp) # 1= M is given as input, 2 = M is an estd matrix (if process_M=0) or a random walk (if process_M=1), 3=Lorenzen, 4=scaled from assessment, 5= M=input M-PAA
data$process_M<-0 # 1= stochastic M for M_model 2 and 3 only, 0= no process error

data$data_simulate <-1 # simulate data, 0=no, 1=yes
data$error_simulate <-1 # simulate process errors, 0=no, 1=yes
data$diet_model <-2 # 1= ddeltadirichlet, 2=ddirichlet
data$flag_nll_diet <- 3 # 1=diet fitted per stomach when age predator unknown , 2=diet fitted per stomach with age predator known, 3=diet fitted to annual average proportions in stomachs, 4=diet fitted to stomach prop averaged over data$xmax years, 5= diet fitted to stomach prop averaged over entire time series 

# For trophic interactions
data$predation_on <- 0 # 0=predation mortality=0, 1=predation on
data$gamma_pref_estim <- 0
data$functional_response <- 1 # Option of functional response, 1=type II, 2=typeIII
data$cons_rate_estim <- 0
data$biomass_other_option <- 1 # Option model biomass_other, 1=cst, 2=surplus production



#### For trophic interactions
if (data$predation_on==1){
  #diet_data<-extract_diet_data_fn("data/diet/diet_stomach_GB_with_empty.csv")
  diet_data<-fake_diet_fn()
  data$age_pred <- diet_data$age_pred
  data$prob_l_given_b <- diet_data$prob_l_given_b # (l,Y,max_B,n_pred)
  data$n_stom <- diet_data$n_stom # (Y,n_pred)
  data$n_stom_max <- max(data$n_stom)
  data$ratio_diet <- diet_data$ratio_diet # (l,Y,sp,n_pred)
  data$length_pred <- diet_data$length_pred

  cons_data<-cons_data_fn("data/diet/diet_consrate_NEUS.csv")
  spring<-cons_data[[1]] # spring survey meansw and bottemp (Y_cons,2,n_pred)
  fall<-cons_data[[2]] # fall survey meansw and bottemp (Y_cons,2,n_pred)
  spring.other<-cons_data[[3]] # spring survey meansw for other food and bottemp (Y_cons,2,n_pred)
  fall.other<-cons_data[[4]] # fall survey meansw for other food and bottemp (Y_cons,2,n_pred)
  data$Y_cons<-cons_data[[5]] # number of years for cons_data


  # Because not the 10 sp yet
  # data$ratio_diet <- array(dim=c(575,data$Y,(data$n_prey),data$n_pred)) #modelled species sp + other food
  # #data$ratio_diet <- diet_data$ratio_diet[,,-(3:10),]
  # for (i in 1:data$n_pred){
  #   data$ratio_diet[,,1:data$sp,i] <- diet_data$ratio_diet[,,1:data$sp,i] # diet cod for prey = cod and had
  #   for (t in 1:data$Y)
  #     data$ratio_diet[,t,(data$sp+1),i] <- apply(diet_data$ratio_diet[,t,(data$sp+1):11,i],1,sum) # diet other species
  # }
  data$spring_cons<-array(dim=c(data$Y_cons,2,data$sp,data$n_pred))
  data$fall_cons<-array(dim=c(data$Y_cons,2,data$sp,data$n_pred))
  data$spring_cons[,,,1]<-spring[,,1:data$sp,1]
  data$fall_cons[,,,1]<-fall[,,1:data$sp,1]
  data$spring_cons_other<-array(dim=c(data$Y_cons,2,data$n_pred))
  data$fall_cons_other<-array(dim=c(data$Y_cons,2,data$n_pred))
  data$spring_cons_other[,,1]<-spring.other[,,1]
  data$fall_cons_other[,,1]<-fall.other[,,1]
  data$flag_cons_rate<-rep(1,data$Y_cons)
  for(t in 1:data$Y_cons){
    for (i in 1:data$sp){
      if (sum(is.na(data$spring_cons[t,,i,1]))>0 || sum(is.na(data$fall_cons[t,,i,1]))>0)
        data$flag_cons_rate[t,i]<-0
    }
  }

} else {
  l=10
  data$age_pred <- array(dim=c(l,data$Y,data$n_pred))
  data$prob_l_given_b <- array(dim=c(l,data$Y,data$max_B,data$n_pred)) # (l,Y,max_B,n_pred)
  data$n_stom <- matrix(2,nrow=data$Y, ncol=data$n_pred) # (Y,n_pred)
  data$n_stom_max <- max(data$n_stom)
  data$ratio_diet <- array(dim=c(l,data$Y,data$sp,data$n_pred)) # (l,Y,sp,n_pred)
  #data$length_pred <- diet_data$length_pred
  data$ratio_diet3=array(1,dim=c(data$Y,data$sp+1,data$max_B,data$n_pred))
  data$Y_cons<-data$Y # number of years for cons_data

  data$spring_cons<-array(dim=c(data$Y_cons,2,data$sp,data$n_pred))
  data$fall_cons<-array(dim=c(data$Y_cons,2,data$sp,data$n_pred))
  data$spring_cons_other<-array(dim=c(data$Y_cons,2,data$n_pred))
  data$fall_cons_other<-array(dim=c(data$Y_cons,2,data$n_pred))
  data$flag_cons_rate<-rep(0,data$Y_cons)
}


if (data$process_F==1){
  for (i in 1:data$sp) {
    if(data$sel_model_F[i]==1){ # need to fix one value so selectivity and random effect not confounded
      init$logit_s_F[data$max_A_catch[i],i]<- 1000 # selectivity of 1 in last catch age
    }
  }
}


#### Prepare options for the map argument of the objective function ####
y_add <- 50
if (data$burnin==1) {
  new_init <- init
  for (i in 1:length(new_init)){
    new_init[[i]] <- burnin.fn(new_init[[i]], names(new_init[i]), y_add)
  }
  new_data <- data
  for (i in 1:length(new_data)){
    new_data[[i]] <- burnin.fn(new_data[[i]], names(new_data[i]), y_add)
  }
  ## special case for t1 and t2
  new_data$t1 <- new_data$t1+y_add
  new_data$t2 <- new_data$t2+y_add
  ## overwrite data and init
  data <- new_data
  init <- new_init
}
source("R/MS_SSM_make.map.fn.r") # file that creates the map argument for objective function
source("R/MS_SSM_make.random.fn.r") # file that creates the random argument for the objective function
map = make.map.fn(data)
random = make.random.fn(data)


# # !!!!! Things that can changed in the map depending on selectivity assumptions
# for (k in 1:data$n_surv_max){
#   for (i in 1:data$sp){
#     if (data$sel_model_surv[i,k]==1){ # !!!!!!!!!!!!!!!!! Check same than in map file
#       for (a in data$max_A_surv[i,k]:data$max_A){
#         init$logit_s_surv[a,i,k]<-1000 #selectivity=1 for the other ages
#       }
#     }
#     if (data$sel_model_F[i]==1){
#       for (a in data$Aplus[i]:data$max_A){ # !!!!!!!!!!!!!!!!! Check same than in map file
#         init$logit_s_F[a,i]<-1000 #selectivity=1 for the other ages
#       }
#     }
#   }
# }
#init$logit_s_surv[3:data$max_A_surv[1,2],1,2]<-1000 #selectivity for cod survey 2
#init$logit_vuln[]=30



#### Create and optimize objective function ####

obj<-MakeADFun(data=data, parameters=init, DLL="MS_SSM", random=random, map=map) 
#save(obj, file="objects_objective_function.RData") #save object of obj for later loading
system.time(opt<-nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 5000, iter.max = 5000,rel.tol=1e-8))) #,sing.tol=1e-20)))

test <- try(
  for(i in seq_len(10)) { # Take 10 extra newton steps
    g <- as.numeric( obj$gr(opt$par) )
    h <- optimHess(opt$par, obj$fn, obj$gr)
    opt$par <- opt$par - solve(h, g)
    opt$objective <- obj$fn(opt$par)
  })

#### Save outputs ####

sd.rep<-sdreport(obj) # save standard errors of estimates
#sd.rep<-sdreport(obj, getReportCovariance = FALSE) # if doesn't need cov matrix
rep<-obj$report(obj$env$last.par.best) # save estimates
# pl <- as.list(sd.rep,"Est")
# plsd <- as.list(sd.rep,"Std")



######################################## Only if model used to simulate data ###################################
#### Simulate data and/or process errors ####

if (data$data_simulate==1 || data$error_simulate==1){
  data2<-obj$simulate(complete=TRUE)
  data2$M_model[]=4
  map = make.map.fn(data2)  
  obj2<-MakeADFun(data=data2, parameters=init, DLL="MS_SSM", random=random, map=map)
  #obj2<-MakeADFun(data=data2, parameters=obj$env$parList(opt$par), DLL="MS_SSM", random=random, map=map)
}
system.time(opt2<-nlminb(obj2$par*0.8, obj2$fn, obj2$gr, control = list(eval.max = 5000, iter.max = 5000,rel.tol=1e-5))) #,sing.tol=1e-20)))

rep2<-obj2$report()
sd.rep2<-sdreport(obj2)


test <- try(
  for(i in seq_len(10)) { # Take 10 extra newton steps
    g <- as.numeric( obj2$gr(opt2$par) )
    h <- optimHess(opt2$par, obj2$fn, obj2$gr)
    opt2$par <- opt2$par - solve(h, g)
    opt2$objective <- obj2$fn(opt2$par)
  })


ch=checkConsistency(obj2,n=600)
summary(ch)
boxplot(t(summary(ch)[[which(names(summary(ch))=="gradientJoint")]]))
abline(h=0,col="red")
boxplot(t(summary(ch)[[which(names(summary(ch))=="gradient")]]))
abline(h=0,col="red")
