extract_diet_data_fn <- function(diet.file){
  
  # From ~\NOAA\Foodhabitsdata\Brian's data\file:///C:/Users/vanessa.trijoulet/Documents/NOAA/Foodhabitsdata/Brian's data/Create_length_to_age_array.R
  # and extract diet
  
  diet<-read.csv(diet.file,header=T)
  diet<-diet[which(is.na(diet$PredGutWeight)==FALSE),] # some NAs in pdgutw
  
  sp_common_names<- c("cod", "haddock", "goosefish", "spinydogfish", "silverhake",
                      "winterskate", "yellowtailflounder", "mackerel", "herring","winterflounder") 
  
  sp_names<-c("GADUS MORHUA","MELANOGRAMMUS AEGLEFINUS","LOPHIUS AMERICANUS","SQUALUS ACANTHIAS","MERLUCCIUS BILINEARIS",
              "LEUCORAJA OCELLATA","LIMANDA FERRUGINEA","SCOMBER SCOMBRUS","CLUPEA HARENGUS","PSEUDOPLEURONECTES AMERICANUS")
  
  prob_l_given_b <- array(0,dim=c(575,data$Y,data$max_B,data$n_pred)) # 575 because only cod for now
  ratio_diet <- array(0,dim=c(575,data$Y,11,data$n_pred))
  n_stom<-array(dim=c(data$Y,data$n_pred))
  length_pred <- array(0,dim=c(575,data$Y,data$n_pred))
  
  for (j in 1:data$n_pred) {
    
    ##### Extract diet by species #####
    
    diet.sp<-diet[diet$predator==sp_names[j],]
    
    name_matrix<-paste("~/NOAA/Foodhabitsdata/Brian's data/GB_",sp_common_names[j],"_length_age_matrix_ageplus.csv",sep="")
    l.a.conv.ageplus<-read.csv(name_matrix,header=FALSE)
    
    
    ### Extract for each year ###
    
    
    for (t in 1:data$Y){
      year<-min(data$year1)+(t-1)
      diet.sp.annual <- diet.sp[diet.sp$Year==year,]
      n_stom[t,j] <- nrow(diet.sp.annual)
      for (n in 1:n_stom[t,j]){
        prob_l_given_b[n,t,,j] <- as.matrix(l.a.conv.ageplus[diet.sp.annual$predLength[n],])
        ratio_diet[n,t,,j] <- unlist(diet.sp.annual[n,8:ncol(diet.sp.annual)])
        length_pred[n,t,j] <- diet.sp.annual$predLength[n]
      }
    }
  }
  return(list(ratio_diet=ratio_diet, prob_l_given_b=prob_l_given_b, n_stom=n_stom, length_pred=length_pred))
}
