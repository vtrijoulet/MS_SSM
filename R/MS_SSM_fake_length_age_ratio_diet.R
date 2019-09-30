fake_diet_fn <- function() {

  sp_common_names<- c("cod", "silverhake", "haddock", "goosefish", "spinydogfish",
                      "winterskate", "yellowtailflounder", "mackerel", "herring","winterflounder") 
  
  n_stom <- array(dim=c(data$Y,data$n_pred))
  n_stom[] = 500
  n_stom_max = max(n_stom)
  fake.ratio_diet <- array(1,dim=c(n_stom_max,data$Y,data$n_prey,data$n_pred)) #just need to be >0 for simulations
  fake.length <- array(0,dim=c(n_stom_max,data$Y,data$n_pred))
  prob_l_given_b <- array(0,dim=c(n_stom_max,data$Y,data$max_B,data$n_pred))
  
  ##### Randomly assign length for ratio_diet ####
  
  
  length=cbind(12:156, 3:76) # min and max for pred from sole = 9:150, modified so starts at age 1 and about same probability for all ages
  # for silver hake 3:76 from age 1
  
  for (j in 1:data$n_pred){
    name_matrix<-paste("data/diet/GB_",sp_common_names[j],"_length_age_matrix_ageplus.csv",sep="")
    l.a.conv.ageplus<-read.csv(name_matrix,header=FALSE)
    for (t in 1:data$Y){
      for (stom in 1:n_stom[t,j]){
        fake.length[stom,t,j]<-sample(length[,j],1)
        prob_l_given_b[stom,t,1:data$Bplus[j],j] <- as.matrix(l.a.conv.ageplus[fake.length[stom,t,j],1:data$Bplus[j]]) 
      }
    }
  }
  
  fake.age <- round(exp(log(fake.length)*lm_coeff[2]+lm_coeff[1]))
  
  # hist(fake.length)
  # hist(fake.age,breaks=0:10)
  
  return(list(ratio_diet=fake.ratio_diet, prob_l_given_b=prob_l_given_b, n_stom=n_stom, length_pred=fake.length, age_pred=fake.age))

}
