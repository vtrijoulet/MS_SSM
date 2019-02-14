### Function that fits the model with shortened data and init from MS_SSM_shorten.fn.R ###


source("R/MS_SSM_shorten.fn.R") # to create new_data and new_init


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

  return(list(data=new_data, init=new_init, map=map, random=random, opt=opt, sd.rep=sd.rep, rep=rep))

}
