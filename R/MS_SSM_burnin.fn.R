library(stringr)

burnin.fn <- function (array, name_array, y_add){
  final_out <- array
  if (is.vector(array)==TRUE){
    if (length(array)==1) {
      if (array==data$Y) final_out<-data$Y+y_add
    } else {
      if (length(array)==data$Y) final_out <- rep(0,data$Y+y_add)
    }
  } else {
    if (sum(dim(array)==data$Y)>0 | sum(dim(array)==data$Y-1)>0){
      dim_array <- which(dim(array)==data$Y)
      dim  <- 1:length(dim(array))
      perm <- c(dim[-dim_array],dim_array) # put years in last dimension
      perm_array <- aperm(array, perm=perm)
      new_dim <- dim(perm_array)
      new_dim[length(new_dim)] <- new_dim[length(new_dim)]+y_add
      mean <- rowMeans(perm_array, dims=length(dim)-1)
      if (str_detect(name_array,"flag")) mean[] <- 0
      out <- array(c(rep(mean, y_add), perm_array), new_dim)
      final_out <- aperm(out, perm=invPerm(perm))
    }
  }
  return(final_out)
}

