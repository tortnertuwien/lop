getLocOut <- function(data, k, alpha=0.5){

  s <- o <- e  <- matrix(NA, nrow=nrow(data), ncol=nrow(data))

  dist <- as.matrix(dist(data))
  core <- knn <- matrix(0, nrow=nrow(data), ncol=nrow(data))

  tmp <- parallel::mclapply(1:nrow(data),
                            FUN=function(x, data, dist, k, alpha){
    return(getLocalDescription(data, center=x, dist, k=k, alpha = alpha))
  }, data=data, dist=dist, k=k, alpha=alpha)

  # lapply(1:nrow(data), FUN=function(x, data, dist, k, alpha){
  #   return(getLocalDescription(data, center=x, dist, k=k, alpha = alpha))
  # }, data=data, dist=dist, k=k, alpha=alpha)

  for(i in 1:nrow(data)){
    e[,i] <- tmp[[i]]$ecd
    s[,i] <- tmp[[i]]$cd
    o[,i] <- tmp[[i]]$od

    core[tmp[[i]]$core,i] <- 1
    knn[tmp[[i]]$knn,i] <- 1
  }

  l <- list(ods = o, cds = s,
            ecds  = e,
            core=core, knn=knn
            )
  class(l) <- "local-outlyingness"
  return( l )
}




