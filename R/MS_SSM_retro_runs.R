## Scripts that runs the retro.fn and plot the results of the retro ##

source("R/MS_SSM_retro.fn.R")


#### Run retro ####

retro_res <- list()
y <- 5 # length of retro run in years

for (n in 1:y){
  retro_res[[paste0(n,"year")]] <- retro.fn(n)
}


#### Make the retro plots ####

# SSB
par(mfrow=c(data$sp,1), oma=c(3,4,1,1), mar=c(0,0,0,0),xpd=NA)
for (i in 1:data$sp){
  plot.time(pred=rep$SSB[,i], se=se$SSB[,i], ylab=paste0(sp.names[i], " SSB"), ylim=c(min(rep$SSB[,i])*0.8, max(rep$SSB[,i])*1.2))
  if(i==data$sp) axis(1)
  for (n in 1:y) lines(retro_res[[n]]$rep$SSB[,i]~years[-((data$Y-n+1):data$Y)], col=n+1)
}


# Recruitment
#par(mfrow=c(data$sp,1), oma=c(3,4,1,1), mar=c(0,0,0,0),xpd=NA)
for (i in 1:data$sp){
  plot.time(pred=rep$NAA[,1,i], se=se$recruits[,i], ylab=paste0(sp.names[i], " recruitment"), ylim=c(min(rep$NAA[,1,i])*0.8, max(rep$NAA[,1,i])*1.2))
  if(i==data$sp) axis(1)
  for (n in 1:y) lines(retro_res[[n]]$rep$NAA[,1,i]~years[-((data$Y-n+1):data$Y)], col=n+1)
}

# mean_Fy
#par(mfrow=c(data$sp,1))
for (i in 1:data$sp){
  plot.time(pred=rep$mean_Fy[,i], se=se$mean_Fy[,i], ylab=paste0(sp.names[i], " mean fishing mortality"), legend=0, ylim=c(min(rep$mean_Fy[,i]*0.6), max(rep$mean_Fy[,i]*1.3)))
  if (i==data$sp) axis(1)
  for (n in 1:y) lines((retro_res[[n]]$rep$mean_Fy[,i])~(years[-((data$Y-n+1):data$Y)]) , col=n+1)
}

# mean_Py
if (data$predation==1){
  for (i in 1:data$sp){
    plot.time(pred=rep$mean_Py[,i], se=se$mean_Py[,i], ylab=paste0("Mean predation mortality on ", sp.names[i]))
    if (i==data$sp) axis(1)
    for (n in 1:y) lines((retro_res[[n]]$rep$mean_Py[,i])~(years[-((data$Y-n+1):data$Y)]) , col=n+1)
  }
}

# Total catch
#par(mfrow=c(data$sp,1), oma=c(3,4,1,1), mar=c(0,0,0,0),xpd=NA)
for (i in 1:data$sp){
  plot.time(pred=rep$aggr_Cw[,i], se=se$aggr_Cw[,i], obs=obs_aggr_Cw[,i], legend=1, ylab=paste0(sp.names[i]," total catch (tons)"),ylim=c(min(rep$aggr_Cw[,i])*0.8, max(rep$aggr_Cw[,i])*1.2))
  if(i==data$sp) axis(1)
  for (n in 1:y) lines((retro_res[[n]]$rep$aggr_Cw[,i])~(years[-((data$Y-n+1):data$Y)]) , col=n+1)
}



