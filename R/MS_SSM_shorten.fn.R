#### Script that cuts down years in a list such as data, init ore rep ####

shorten.fn <- function(list, x=NULL, range=NULL) { # x=years (e.g. 5) to remove final years (e.g. retro), range=range (1:50) of years to remove for any other shortening
  
  if (missing(x) & missing(range)) stop("Need to specified x or range")
  if (!missing(x) & !missing(range)) stop("Cannot specified both x and range")
  
  ## Find dimension with year index
  dim_list<- list()
  for (i in 1: length(list)){
    if (is.vector(list[[i]])==TRUE){
      if (length(list[[i]])==1) {
        if (is.nan(list[[i]])) dim_list[i] <- 0 else
          if (list[[i]]==data$Y | list[[i]]==data$Y-1) dim_list[i] <- 1000 else dim_list[i] <- 0
      } else {
        if (length(list[[i]])==data$Y | length(list[[i]])==data$Y-1) dim_list[i] <- 2000 else dim_list[i] <- 0
      }
    } else {
      if (sum(dim(list[[i]])==data$Y)>0) dim_list[i] <- which(dim(list[[i]])==data$Y) else 
        if (sum(dim(list[[i]])==data$Y-1)>0) dim_list[i] <- which(dim(list[[i]])==data$Y-1) else dim_list[i] <- 0
    }
  }
  
  
  ## Create new init and data with less years
  new_list <- list
  for (i in 1:length(new_list)){
    if (!missing(x)) { # years to remove from the end (retro)
      if (dim_list[[i]]==1000) new_list[[i]]<-data$Y-x
      if (dim_list[[i]]==2000) new_list[[i]]<-new_list[[i]][-((length(new_list[[i]])-x+1):length(new_list[[i]]))]
      if (dim_list[[i]]==1){
        dim <- length(dim(new_list[[i]]))
        comma <-paste((rep(",",dim-1)), collapse ="")
        expr <- paste0("new_list[[i]][-((nrow(new_list[[i]])-x+1):nrow(new_list[[i]]))",comma,",drop=FALSE]")
        new_list[[i]] <- eval(parse(text=expr))
      }
      if (dim_list[[i]]==2){
        dim <- length(dim(new_list[[i]]))
        comma <-paste((rep(",",dim-2)), collapse ="")
        expr <- paste0("new_list[[i]][,-((nrow(new_list[[i]])-x+1):nrow(new_list[[i]]))",comma,",drop=FALSE]")
        new_list[[i]] <- eval(parse(text=expr))
      }
    } else { # age range instead (e.g. remove burn-in years)
      if (dim_list[[i]]==1000) new_list[[i]]<-data$Y-length(range)
      if (dim_list[[i]]==2000) new_list[[i]]<-new_list[[i]][-range]
      if (dim_list[[i]]==1){
        dim <- length(dim(new_list[[i]]))
        comma <-paste((rep(",",dim-1)), collapse ="")
        expr <- paste0("new_list[[i]][-range",comma,",drop=FALSE]")
        new_list[[i]] <- eval(parse(text=expr))
      }
      if (dim_list[[i]]==2){
        dim <- length(dim(new_list[[i]]))
        comma <-paste((rep(",",dim-2)), collapse ="")
        expr <- paste0("new_list[[i]][,-range",comma,",drop=FALSE]")
        new_list[[i]] <- eval(parse(text=expr))
      }
      
    }
  }
  return(new_list)
}
