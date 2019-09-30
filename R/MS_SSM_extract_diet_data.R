extract_diet_data_fn <- function(diet.file){
  
  # From ~\Foodhabitsdata\Brian's data\Create_length_to_age_array.R
  # and extract diet
  
  diet<-read.csv(diet.file,header=T)
  diet<-diet[which(is.na(diet$PredGutWeight)==FALSE),] # some NAs in pdgutw
  
  sp_common_names<- c("cod", "silverhake", "herring", "haddock", "goosefish", 
                      "spinydogfish","winterskate", "yellowtailflounder", "mackerel","winterflounder") 
  
  sp_names<-c("GADUS MORHUA","MERLUCCIUS BILINEARIS","CLUPEA HARENGUS","MELANOGRAMMUS AEGLEFINUS","LOPHIUS AMERICANUS",
              "SQUALUS ACANTHIAS","LEUCORAJA OCELLATA","LIMANDA FERRUGINEA","SCOMBER SCOMBRUS","PSEUDOPLEURONECTES AMERICANUS")
  
  ## Estimate max number of annual stomachs for predator of interest
  n_stom_max=c()
  for (j in 1:length(sp_names)) {
    diet.sp<-diet[diet$predator==sp_names[j],]
    n_stom=c()
    for (t in 1:data$Y){
      year<-min(data$year1)+(t-1)
      diet.sp.annual <- diet.sp[diet.sp$Year==year,]
      n_stom <- c(n_stom,nrow(diet.sp.annual))
    }
    n_stom_max <- c(n_stom_max,max(n_stom))
  }
  max_n_stom <- max(n_stom_max[1:data$n_pred])
  
  ## Create object with max_n_stom
  prob_l_given_b <- array(0,dim=c(max_n_stom,data$Y,data$max_B,data$n_pred)) 
  ratio_diet <- array(0,dim=c(max_n_stom,data$Y,length(sp_common_names)+1,data$n_pred))
  n_stom<-array(dim=c(data$Y,data$n_pred))
  length_pred <- array(0,dim=c(max_n_stom,data$Y,data$n_pred))
  
  
  
  for (j in 1:data$n_pred) {
    
    ##### Extract diet by species #####
    
    diet.sp<-diet[diet$predator==sp_names[j],]
    
    name_matrix<-paste("data/diet/GB_",sp_common_names[j],"_length_age_matrix_ageplus.csv",sep="")
    l.a.conv.ageplus<-read.csv(name_matrix,header=FALSE)
    
    
    ### Extract for each year ###
    
    
    for (t in 1:data$Y){
      year<-min(data$year1)+(t-1)
      diet.sp.annual <- diet.sp[diet.sp$Year==year,]
      n_stom[t,j] <- nrow(diet.sp.annual)
      for (n in 1:n_stom[t,j]){
        prob_l_given_b[n,t,1:data$Aplus[j],j] <- as.matrix(l.a.conv.ageplus[diet.sp.annual$predLength[n],])
        ratio_diet[n,t,,j] <- unlist(diet.sp.annual[n,8:ncol(diet.sp.annual)])
        length_pred[n,t,j] <- diet.sp.annual$predLength[n]
      }
    }
  }
  return(list(ratio_diet=ratio_diet, prob_l_given_b=prob_l_given_b, n_stom=n_stom, length_pred=length_pred))
}
