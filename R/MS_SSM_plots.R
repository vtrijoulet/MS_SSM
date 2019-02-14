years<-min(data$year1):data$lastyear
Y1<-min(data$year1)
sp.names<-attr(data,"sp_names")
z.stat <- 1.96

name.folder<-paste(paste(sp.names,collapse ="_"),"_",Y1,"-",data$lastyear,"_rec",data$recruit_model,
                "_M",data$M_model,"_err-rec",data$process_rec,"_err-survi",data$process_survival,
                "_err-M",data$process_M,"_pred",data$predation_on,"_Type",data$functional_response+1,
                "_sizepref",data$gamma_pref_estim,"_consrate",data$cons_rate_estim,"_flagdiet",data$flag_nll_diet,sep="")[1]

mypath<-getwd()
new.path<-paste(mypath,"Results",name.folder,sep="/")
dir.create(new.path,recursive=TRUE)
setwd(new.path)

#### Calculate CI ####
names.rep.all <- names(sd.rep$value)
names.rep <- unique(names.rep.all)
se <- list()
k=0
for (i in 1:length(names.rep)){
  k=k+1
  dim <- dim(rep[[names.rep[i]]])
  idx <- as.vector(array(k:(k+prod(dim)-1),dim=dim))
  se[[names.rep[i]]] <- array(sd.rep$sd[idx], dim=dim)
  k=k+prod(dim)-1
}


# #### Calculate mean mortality rates and corresponding CI ####
# 
# ## Function to estimate se for mean of adereported objects
# require(numDeriv)
# se.ADrep <- function (object, sp, FUN){
#   phi <- function(x) FUN(x)
#   ix <- obj$env$ADreportIndex()[[object]][,,sp]
#   covx <- sd.rep$cov[ix, ix]
#   x <- rep[[object]][,,sp]
#   J <- jacobian(phi, x) ## Derivative of phi at x
#   covPhix <- J %*% covx %*% t(J) ## Covariance of phi(x)
#   return(sqrt(diag(covPhix))) # se around phi(x)
# }
# 
# ## F
# mean.F<-array(dim=c(data$Y,data$sp))
# se.mean.F<-array(0,dim=c(data$Y,data$sp))
# mean.FAA<-array(dim=c(data$max_A,data$sp))
# se.mean.FAA<-array(0,dim=c(data$max_A,data$sp))
# for (i in 1:data$sp){
#   mean.F[,i]<-apply(rep$F[,data$min_A_catch[i]:data$max_A_catch[i],i],1,mean)
#   se.mean.F[,i] <- se.ADrep("F",i,rowMeans)
#   mean.FAA[data$min_A_catch[i]:data$max_A_catch[i],i]<-apply(rep$F[,data$min_A_catch[i]:data$max_A_catch[i],i],2,mean)
#   se.mean.FAA[,i] <- se.ADrep("F",i,colMeans)
# }
# 
# ## M 
# mean.MAA<-array(dim=c(data$max_A,data$sp))
# se.mean.MAA<-array(0,dim=c(data$max_A,data$sp))
# for (i in 1:data$sp){
#   mean.MAA[1:data$Aplus[i],i]<-apply(rep$MAA[,1:data$Aplus[i],i],2,mean)
#   se.mean.MAA[,i] <- se.ADrep("MAA",i,colMeans)
# }
# 
# ## P
# if (data$predation_on==1){
#   mean.PAA<-array(dim=c(data$max_A,data$sp))
#   se.mean.PAA<-array(0,dim=c(data$max_A,data$sp))
#   mean.PAA.y<-array(dim=c(data$Y,data$sp))
#   se.mean.PAA.y<-array(dim=c(data$Y,data$sp))
#   for (i in 1:data$sp){
#     mean.PAA[,i]<-apply(rep$PAA[,,i],2,mean)
#     se.mean.PAA[,i]<-se.ADrep("PAA",i,colMeans)
#     mean.PAA.y[,i]<-apply(rep$PAA[,1:data$Aplus[i],i],1,mean)
#     se.mean.PAA.y[,i]<-se.ADrep("PAA",i,rowMeans)
#   }
# }
# 
# ## Z
# mean.Z<-array(dim=c(data$max_A,data$sp))
# se.mean.Z<-array(dim=c(data$max_A,data$sp))
# mean.Z.y<-array(dim=c(data$Y,data$sp))
# se.mean.Z.y<-array(dim=c(data$Y,data$sp))
# for (i in 1:data$sp){
#   mean.Z[,i]<-apply(rep$Z[,,i],2,mean)
#   se.mean.Z[,i]<-se.ADrep("Z",i,colMeans)
#   mean.Z.y[,i]<-apply(rep$Z[,1:data$Aplus[i],i],1,mean)
#   se.mean.Z.y[,i]<-se.ADrep("Z",i,rowMeans)
# }


#### Plot functions ####

plot.time <- function(pred, obs=NULL, se=0, legend=0, xaxt="n", xlab="", ylab="", main="", ylim=c(min(na.omit(c((pred-z.stat*se),obs))),max(na.omit(c((pred+z.stat*se),obs)))) ){
  plot(pred~years, ylim=ylim, ylab=ylab, type="l", main=main, xlab=xlab, xaxt=xaxt)
  if (!is.null(obs)) points(obs~years, pch=16, col="red")
  if (!is.null(se)) polygon(x=c(years,rev(years)), y=c(pred-z.stat*se,rev(pred+z.stat*se)), col=rgb(0,0,0,alpha=0.3), border = NA)
  if (legend==1 & i==1) legend("topleft",legend=c("Observed","Predicted"),pch=c(16,NA),lty=c(NA,1),col=c("red", "black"),bty="n", cex=0.8)
}
plot.age <- function(pred, age, se=0, legend=0, xaxt="n", xlab="", ylab="", main="", ylim=c(min(na.omit((pred-z.stat*se))),max(na.omit((pred+z.stat*se)))) ){
  plot(y=pred, x=age, ylim=ylim, ylab=ylab, type="l", xaxt=xaxt, main=main, xlab=xlab, xlim=c(1, data$max_A))
  if (!is.null(se)) polygon(x=c(age,rev(age)), y=c(pred-z.stat*se,rev(pred+z.stat*se)), col=rgb(0,0,0,alpha=0.3), border = NA)
  #if (legend==1 & i==1) legend("topright",legend=c("Observed","Predicted"),pch=c(16,NA),lty=c(NA,1),col=c("red", "black"),bty="n")
}


#### Start pdf output file ####
pdf(file="Figures_fit.pdf")

obs_aggr_Cw<-data$obs_aggr_Cw
for (i in 1:data$sp)
  for (t in 1:data$Y)
    if (data$obs_aggr_Cw[t,i]==-999)
      obs_aggr_Cw[t,i]<-NA

par(mfrow=c(data$sp,1), oma=c(3,4,1,1), mar=c(0,0,0,0),xpd=NA)
for (i in 1:data$sp){
  plot.time(pred=rep$aggr_Cw[,i], obs=obs_aggr_Cw[,i], se=se$aggr_Cw[,i], legend=1, ylab=paste0(sp.names[i]," total catch (tons)"))
  if(i==data$sp) axis(1)
}



obs_aggr_I<-data$obs_aggr_I
for (j in 1:data$n_surv[i])
  for (i in 1:data$sp)
    for (t in 1:data$Y)
      if (data$obs_aggr_I[t,i,j]==-999)
        obs_aggr_I[t,i,j]<-NA

par(mfrow=c(data$sp,data$n_surv_max), oma=c(3,4,2,1), mar=c(0,1,0,1),xpd=NA)
if (data$n_surv_max>4) {par(mfrow=c(data$sp,3), oma=c(3,4,2,1), mar=c(0,1,0,1),xpd=NA)
  } else {par(mfrow=c(data$sp,data$n_surv_max), oma=c(3,4,2,1), mar=c(0,1,0,1),xpd=NA)}
for (i in 1:data$sp){
  for (j in 1:data$n_surv[i]){
    plot.time(pred=rep$aggr_I[,i,j], obs=obs_aggr_I[,i,j], se=se$aggr_I[,i,j], legend=1)
    if (i==1) mtext(paste0("Survey ",j),side=3, line=0.5, cex=0.7)
    if(i==data$sp) axis(1)
    if(j==1) mtext(paste0(sp.names[i]," total abundance indices"), side=2, line=3, cex=0.6)
  }
}


for (j in 1:data$sp){
  if (length(data$min_A_catch[j]:data$max_A_catch[j])>4)
    par(mfrow=c(ceiling(length(data$min_A_catch[j]:data$max_A_catch[j])/3),3), oma=c(2,2,1,1), mar=c(2,2,2,2),xpd=NA)
  else 
    par(mfrow=c(ceiling(length(data$min_A_catch[j]:data$max_A_catch[j])/2),2), oma=c(2,2,1,1), mar=c(2,2,2,2),xpd=NA)
  for (i in data$min_A_catch[j]:data$max_A_catch[j]){
    plot.time(pred=rep$prop_Caa[,i,j], obs=data$obs_prop_Caa[,i,j], xaxt="s", legend=1, main=paste0("Age ",i), ylab=paste(sp.names[j],"age comp catch",sep=" "))
  }
  #for (n in 1:(12-length(data$min_A_catch[j]:data$max_A_catch[j])))
    #plot.new()
}


#par(mfrow=c(4,3), oma=c(2,2,1,1), mar=c(2,2,2,2),xpd=NA)
for (j in 1:data$sp){
  for (k in 1:data$n_surv[j]){
    if (length(data$min_A_surv[j,k]:data$max_A_surv[j,k])>4)
      par(mfrow=c(ceiling(length(data$min_A_surv[j,k]:data$max_A_surv[j,k])/3),3), oma=c(2,2,1,1), mar=c(2,2,2,2),xpd=NA)
    else 
      par(mfrow=c(ceiling(length(data$min_A_surv[j,k]:data$max_A_surv[j,k])/2),2), oma=c(2,2,1,1), mar=c(2,2,2,2),xpd=NA)
    if (length(data$min_A_surv[j,k]:data$max_A_surv[j,k])==1) par(mfrow=c(1,1), oma=c(2,2,1,1), mar=c(2,2,2,2),xpd=NA)
    for (i in data$min_A_surv[j,k]:data$max_A_surv[j,k]){
      plot.time(pred=rep$prop_Iaa[,i,j,k], obs=data$obs_prop_Iaa[,i,j,k], xaxt="s", legend=1, main=paste0("Age ",i), ylab=paste(sp.names[j]," age comp survey",k,sep=" "))
    }
    #for (n in 1:(12-length(data$min_A_surv[j,k]:data$max_A_surv[j,k])))
      #plot.new()
  }
}


# Plot SSB
par(mfrow=c(data$sp,1), oma=c(3,4,1,1), mar=c(0,0,0,0),xpd=NA)
for (i in 1:data$sp){
  plot.time(pred=rep$SSB[,i], se=se$SSB[,i], ylab=paste0(sp.names[i], " SSB"))
  if(i==data$sp) axis(1)
}


#plot recruits
#par(mfrow=c(data$sp,1), oma=c(3,4,1,1), mar=c(0,0,0,0),xpd=NA)
for (i in 1:data$sp){
  plot.time(pred=rep$NAA[,1,i], se=se$recruits[,i], ylab=paste0(sp.names[i], " recruitment"))
  if(i==data$sp) axis(1)
}


par(mfrow=c(data$sp,1), oma=c(2,2,1,1), mar=c(2,2,2,2),xpd=NA)
for (i in 1:data$sp){
  plot(rep$NAA[2:data$Y,1,i]~rep$SSB[1:(data$Y-1),i], type="l", xlab=paste0(sp.names[i]," SSB"), ylab=paste0(sp.names[i]," recruitment"))
  text(y=rep$NAA[2:data$Y,1,i], x=rep$SSB[1:(data$Y-1),i], label=years[-data$Y], col="red")
}


# Plot s_surv
par(mfrow=c(data$sp,data$n_surv_max), oma=c(3,4,2,1), mar=c(0,1,0,1),xpd=NA)
for (i in 1:data$sp){
  for (j in 1:data$n_surv[i]){
    plot.age(pred=rep$s_surv[data$min_A_surv[i,j]:data$max_A_surv[i,j],i,j], age=data$min_A_surv[i,j]:data$max_A_surv[i,j])
    if (i==1) mtext(paste0("Survey ",j),side=3, line=0.5, cex=0.7)
    if (j==1) mtext(paste0("Survey selectivity on ", sp.names[i]), side=2, line=3, cex=0.6)
  if (i==data$sp) axis(1)
  }
}

# Plot s_F
par(mfrow=c(data$sp,1), oma=c(3,4,1,1), mar=c(0,0,0,0),xpd=NA)
for (i in 1:data$sp){
  plot.age(pred=rep$s_F[data$min_A_catch[i]:data$max_A_catch[i],i], age=data$min_A_catch[i]:data$max_A_catch[i], ylab=paste0("Fishing selectivity on ", sp.names[i]))
  if (i==data$sp) axis(1)
}


# Plot F
#par(mfrow=c(data$sp,1))
for (i in 1:data$sp){
  plot.time(pred=rep$mean_Fy[,i], se=se$mean_Fy[,i], ylab=paste0(sp.names[i], " mean fishing mortality"), legend=0)
  if (i==data$sp) axis(1)
}

#par(mfrow=c(data$sp,1))
for (i in 1:data$sp){
  plot.age(pred=rep$mean_FAA[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i], se=se$mean_FAA[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i], age=data$min_Fbar_idx[i]:data$max_Fbar_idx[i], ylab=paste0("Mean fishing mortality"))
  if (i==data$sp) axis(1)
}

#par(mfrow=c(data$sp,1))
for (i in 1:data$sp) {
  plot.age(pred=rep$F[1,1:data$Aplus[i],i], age=1:data$Aplus[i], ylab="Fishing mortality in all years", ylim=c(0,max(rep$F[,,i])))
  for(t in 2:data$Y){
    lines(rep$F[t,1:data$Aplus[i],i],col=t)
  }
  if (i==data$sp) axis(1)
}

#par(mfrow=c(data$sp,1))
for (i in 1:data$sp) {
  plot.time(pred=rep$F[,1,i], ylim=c(min(rep$F[,,i]),max(rep$F[,,i])), ylab=paste0("Fishing mortality at all ages on ", sp.names[i]))
  for(a in 2:data$Aplus[i]){
    lines(rep$F[,a,i]~years,col=a)
  }
  if (i==1) legend("topleft", lty=1, col=1:data$max_A, legend=1:data$max_A, bty="n", cex=0.7, ncol=2)
  if (i==data$sp) axis(1)
}



#plot M
#par(mfrow=c(data$sp,1))
for (i in 1:data$sp) {
  plot.age(pred=rep$MAA[1,1:data$Aplus[i],i], age=1:data$Aplus[i], ylab="Natural mortality in all years", ylim=c(min(rep$MAA),max(rep$MAA)))
  for(t in 2:data$Y){
    lines(rep$MAA[t,1:data$Aplus[i],i],col=t)
  }
  if (i==data$sp) axis(1)
}

#par(mfrow=c(data$sp,1))
for (i in 1:data$sp){
  plot.time(pred=rep$mean_My[,i], se=se$mean_My[,i], ylab=paste0(sp.names[i], " mean natural mortality"), legend=0)
  if (i==data$sp) axis(1)
}

#par(mfrow=c(data$sp,1))
for (i in 1:data$sp){
  plot.age(pred=rep$mean_MAA[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i], se=se$mean_MAA[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i], age=data$min_Fbar_idx[i]:data$max_Fbar_idx[i], ylab="Mean natural mortality")
}



#plot PAA
#par(mfrow=c(data$sp,1))
if (data$predation_on==1){
  for (i in 1:data$sp){
    plot.age(pred=rep$mean_PAA[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i], se=se$mean_PAA[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i], age=data$min_Fbar_idx[i]:data$max_Fbar_idx[i], ylab="Mean predation mortality")
    if (i==data$sp) axis(1)
  }
  for (i in 1:data$sp){
    plot.time(pred=rep$mean_Py[,i], se=se$mean_Py[,i], ylab=paste0("Mean predation mortality on ", sp.names[i]))
    if (i==data$sp) axis(1)
  }
  
  #par(mfrow=c(data$sp,1))
  for (i in 1:data$sp) {
    plot.age(pred=rep$PAA[1,1:data$Aplus[i],i], age=1:data$Aplus[i], ylab="Predation mortality in all years", ylim=c(0,max(rep$PAA[,,i])))
    for(t in 2:data$Y){
      lines(rep$PAA[t,1:data$Aplus[i],i],col=t)
    }
    if (i==data$sp) axis(1)
  }
  
  #par(mfrow=c(data$sp,1))
  for (i in 1:data$sp) {
    plot.time(pred=rep$PAA[,1,i], ylim=c(min(rep$PAA[,,i]),max(rep$PAA[,,i])), ylab=paste0("Predation mortality at all ages on ", sp.names[i]))
    for(a in 2:data$Aplus[i]){
      lines(rep$PAA[,a,i]~years,col=a)
    }
    if (i==1) legend("topleft", lty=1, col=1:data$max_A, legend=1:data$max_A, bty="n", cex=0.7, ncol=2)
    if (i==data$sp) axis(1)
  }
}


#plot Z

#par(mfrow=c(data$sp,1))
for (i in 1:data$sp){
  plot.age(pred=rep$mean_ZAA[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i], se=se$mean_ZAA[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i], age=data$min_Fbar_idx[i]:data$max_Fbar_idx[i], ylab="Mean total mortality")
  if (i==data$sp) axis(1)
}
for (i in 1:data$sp){
  plot.time(pred=rep$mean_Zy[,i], se=se$mean_Zy[,i], ylab=paste0("Mean total mortality on ", sp.names[i]))
  if (i==data$sp) axis(1)
}
#par(mfrow=c(data$sp,1))
for (i in 1:data$sp) {
  plot.age(pred=rep$Z[1,1:data$Aplus[i],i], age=1:data$Aplus[i], ylab="Total mortality in all years", ylim=c(min(rep$Z),max(rep$Z)))
  for(t in 2:data$Y){
    lines(rep$Z[t,1:data$Aplus[i],i],col=t)
  }
  if (i==data$sp) axis(1)
}


#Area plot mortalities
col<-c("deepskyblue3","midnightblue","darkolivegreen3")
if (data$predation_on==0) {
  rep$mean_PAA=rep$mean_MAA
  rep$mean_PAA[]=0
}
M=rep$mean_PAA+rep$mean_MAA
#par(mfrow=c(data$sp,1))
for (i in 1:data$sp){
  age=(data$min_Fbar_idx[i]:data$max_Fbar_idx[i])
  plot(rep$mean_ZAA[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i]~age,xaxt="n", type="n", xlab="",ylab=paste0("Total mortality on ", sp.names[i]),ylim=c(0,max(rep$mean_ZAA[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i])))
  polygon(x=c(data$min_Fbar_idx[i]:data$max_Fbar_idx[i],rev(data$min_Fbar_idx[i]:data$max_Fbar_idx[i])), c(rep$mean_PAA[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i], rev(rep(0,length(data$min_Fbar_idx[i]:data$max_Fbar_idx[i])))),col=col[1],border=NA)
  polygon(x=c(data$min_Fbar_idx[i]:data$max_Fbar_idx[i],rev(data$min_Fbar_idx[i]:data$max_Fbar_idx[i])), c(M[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i], rev(rep$mean_PAA[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i])),col=col[2],border=NA)
  polygon(x=c(data$min_Fbar_idx[i]:data$max_Fbar_idx[i],rev(data$min_Fbar_idx[i]:data$max_Fbar_idx[i])), c(rep$mean_ZAA[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i], rev(M[data$min_Fbar_idx[i]:data$max_Fbar_idx[i],i])),col=col[3],border=NA)
  if (i==1) if (data$predation==1) legend("topleft", legend=rev(c("P", "M", "F")), fill=rev(col), border=NA, bty="n") else legend("topleft", legend=rev(c("M", "F")), fill=rev(col), border=NA, bty="n")
  if (i==data$sp) axis(1)
}

if (data$predation_on==0) {
  rep$mean_Py=rep$mean_My
  rep$mean_Py[]=0
}
My=rep$mean_Py+rep$mean_My
#par(mfrow=c(data$sp,1))
for (i in 1:data$sp){
  plot(rep$mean_Zy[,i]~years,xaxt="n", type="n", xlab="",ylab=paste0("Total mortality on ", sp.names[i]),ylim=c(0,max(rep$mean_Zy[,i])))
  polygon(x=c(years,rev(years)), c(rep$mean_Py[,i], rev(rep(0,data$Y))),col=col[1],border=NA)
  polygon(x=c(years,rev(years)), c(My[,i], rev(rep$mean_Py[,i])),col=col[2],border=NA)
  polygon(x=c(years,rev(years)), c(rep$mean_Zy[,i], rev(My[,i])),col=col[3],border=NA)
  if (i==1) if (data$predation==1) legend("topleft", legend=rev(c("P", "M", "F")), fill=rev(col), border=NA, bty="n") else legend("topleft", legend=rev(c("M", "F")), fill=rev(col), border=NA, bty="n")
  if (i==data$sp) axis(1)
}



# # Plot cons_rate
# if (data$cons_rate_estim==1){
#   par(mfrow=c(data$sp,data$n_pred))
#   for (j in 1:data$n_pred)
#     for (i in 1:data$sp){
#       plot(rep$cons_rate[i,,j],main=sp.names[i],ylab="Consumption rate (kg/predator/year)",xlab=paste("Age",sp.names[j], "predator",sep=" "),type="l")
#     }
# }


dev.off() ## end of pdf file


#### Save estimates #####
# outputs<-cbind(sd.rep$value,sd.rep$sd)
# colnames(outputs)<-c("Estimate","Std.Error")
# fixed.param<-cbind(sd.rep$par.fixed)
# #random.param<-cbind(sd.rep$par.random)
# write.csv(outputs,"sd_rep.csv")
# write.csv(fixed.param,"fixed_par.csv")
# write.csv(rep,"rep.csv")
# sink("outtext.txt") 
# lapply(rep, print) 
# sink() 
save(data,file="data.RData")
save(init,file="init.RData")
save(rep,file="rep.RData")
save(sd.rep,file="sd.rep.RData")
save(map,file="map.RData")
save(opt,file="opt.RData")
#save.image("MS_SSM_workspace.RData")