### Function that fits the model with shortened data and init from MS_SSM_shorten.fn.R ###


source("R/MS_SSM_shorten.fn.R") # to create new_data and new_init
library(numDeriv)


retro.fn <- function(n){
  
  new_data <- shorten.fn(data, x=n) # creation of new_data 
  new_init <- shorten.fn(init,x=n) # creation of new_init
  map = make.map.fn(new_data)
  random = make.random.fn(new_data)
  
  obj<-MakeADFun(data=new_data, parameters=new_init, DLL="MS_SSM", random=random, map=map) 
  system.time(opt<-nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 5000, iter.max = 5000,rel.tol=1e-8))) #,sing.tol=1e-20)))
  
  test <- try(
    for(i in seq_len(10)) { # Take 10 extra newton steps
      g <- as.numeric( obj$gr(opt$par) )
      h <- optimHess(opt$par, obj$fn, obj$gr)
      opt$par <- opt$par - solve(h, g)
      opt$objective <- obj$fn(opt$par)
    })
  
  sd.rep<-sdreport(obj) # save standard errors of estimates
  rep<-obj$report(obj$env$last.par.best) 
  
  se.ADrep <- function (object, sp, FUN){
    phi <- function(x) FUN(x)
    ix <- obj$env$ADreportIndex()[[object]][,,sp]
    covx <- sd.rep$cov[ix, ix]
    x <- rep[[object]][,,sp]
    J <- jacobian(phi, x) ## Derivative of phi at x
    covPhix <- J %*% covx %*% t(J) ## Covariance of phi(x)
    return(sqrt(diag(covPhix))) # se around phi(x)
  }
  mean.F<-array(dim=c(new_data$Y,new_data$sp))
  se.mean.F<-array(0,dim=c(new_data$Y,new_data$sp))
  for (i in 1:new_data$sp){
    mean.F[,i]<-apply(rep$F[,new_data$min_A_catch[i]:new_data$max_A_catch[i],i],1,mean)
    se.mean.F[,i] <- se.ADrep("F",i,rowMeans)
  }
  if (new_data$predation_on==1){
    mean.PAA.y<-array(dim=c(new_data$Y,new_data$sp))
    se.mean.PAA.y<-array(dim=c(new_data$Y,new_data$sp))
    for (i in 1:new_data$sp){
      mean.PAA.y[,i]<-apply(rep$PAA[,1:new_data$Aplus[i],i],1,mean)
      se.mean.PAA.y[,i]<-se.ADrep("PAA",i,rowMeans)
    }
  }
  if (new_data$predation_on==1){
    return(list(data=new_data, init=new_init, map=map, random=random, opt=opt, sd.rep=sd.rep, rep=rep, mean.F=mean.F, se.mean.F=se.mean.F, mean.PAA.y=mean.PAA.y, se.mean.PAA.y=se.mean.PAA.y))
  } else {
    return(list(data=new_data, init=new_init, map=map, random=random, opt=opt, sd.rep=sd.rep, rep=rep, mean.F=mean.F, se.mean.F=se.mean.F))
  }
}
