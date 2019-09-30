##### Prepare data for observed consumption rates #####

cons_data_fn <- function (data.file.diet2){
  
  data<-read.csv(data.file.diet2)
  
  # cf.("Consumption rates/Make meansw and bottemp as data_per prey species.R")
  
  
  # Separate other food from consumption for 10 modelled sp
  data2<-data[which(data$prey!="OTHER"),]
  data.other<-data[which(data$prey=="OTHER"),]
  
  # ignore data when <5 tows
  data3 <- data2
  data.other3 <- data.other
  for (n in 1:nrow(data3)){
    if (data3$num_tows[n]<5){
      data3$meansw[n]<- NA
    }
  }  
  for (n in 1:nrow(data.other3)){
    if (data.other3$num_tows[n]<5 || is.na(data.other3$num_tows[n])==TRUE){
      data.other3$meansw[n]<- NA
    }
  }
  
  # ignore data when <20 stomachs
  data4 <- data3
  data.other4 <- data.other3
  for (n in 1:nrow(data4)){
    if (data4$nstom[n]<20){
      data4$meansw[n]<- NA
    }
  }
  for (n in 1:nrow(data.other4)){
    if (data.other4$nstom[n]<20){
      data.other4$meansw[n]<- NA
    }
  }
  
  #### Make data so rows=years and columns = meansw (summed over species) and bottemp for each survey
  years<-sort(unique(data4$year))
  pred<-unique(data4$svspp)
  Y<-length(years)
  # sp_common_names<- c("cod", "haddock", "goosefish", "spinydogfish", "silverhake",
  #                     "winterskate", "yellowtailflounder", "mackerel", "herring","winterflounder") 
  # sp<-c("GADUS MORHUA","MELANOGRAMMUS AEGLEFINUS","LOPHIUS AMERICANUS","SQUALUS ACANTHIAS","MERLUCCIUS BILINEARIS",
  #       "LEUCORAJA OCELLATA","LIMANDA FERRUGINEA","SCOMBER SCOMBRUS","CLUPEA HARENGUS","PSEUDOPLEURONECTES AMERICANUS")
  svspp_order<-c(73,74,197,15,72,
                 23,105,121,32,106)
  prey<-c("GADMOR","MELAEG","LOPAME","SQUACA","MERBIL","LIMFER","SCOSCO","CLUHAR","PLEAME")
  n_pred<-length(svspp_order)
  n_prey<-length(prey)
  spring<-array(0,dim=c(Y,2,length(prey),n_pred))
  fall<-array(0,dim=c(Y,2,length(prey),n_pred))
  spring.other<-array(0,dim=c(Y,2,n_pred))
  fall.other<-array(0,dim=c(Y,2,n_pred))
  
  for (j in 1:n_pred){
    species<-data4[which(data4$svspp==svspp_order[j]),]
    species.other<-data.other4[which(data.other4$svspp==svspp_order[j]),]
    for (t in 1:Y){
      annual<-species[which(species$year==years[t]),]
      if (length(annual$meansw[which(annual$seacat=="SPRING")])==0){
        spring[t,,,j]<-NA
      } else {
        for (i in 1:n_prey){
          annual_prey<-annual[which(annual$prey==prey[i]),]
          spring[t,1,i,j]<-annual_prey$meansw[which(annual_prey$seacat=="SPRING")]
          spring[t,2,i,j]<-annual_prey$bottemp[which(annual_prey$seacat=="SPRING")]
        }
      }
      if (length(annual$meansw[which(annual$seacat=="FALL")])==0){
        fall[t,,,j]<-NA
      } else {
        for (i in 1:n_prey){
          annual_prey<-annual[which(annual$prey==prey[i]),]
          fall[t,1,i,j]<-annual_prey$meansw[which(annual_prey$seacat=="FALL")]
          fall[t,2,i,j]<-annual_prey$bottemp[which(annual_prey$seacat=="FALL")]
        }
      }
      annual.other<-species.other[which(species.other$year==years[t]),]
      if (length(annual.other$meansw[which(annual.other$seacat=="SPRING")])==0){
        spring.other[t,,j]<-NA
      } else {
        spring.other[t,1,j]<-annual.other$meansw[which(annual.other$seacat=="SPRING")]
        spring.other[t,2,j]<-annual.other$bottemp[which(annual.other$seacat=="SPRING")]
      }
      if (length(annual.other$meansw[which(annual.other$seacat=="FALL")])==0){
        fall.other[t,,j]<-NA
      } else{
        fall.other[t,1,j]<-annual.other$meansw[which(annual.other$seacat=="FALL")]
        fall.other[t,2,j]<-annual.other$bottemp[which(annual.other$seacat=="FALL")]
      }
    }
  }
  
  
  
  return(list(spring,fall,spring.other,fall.other,Y))
  
}