evaluateCores <- function(cores, groups, ...){

  ug <- levels(groups)
  qm <- extractQualityMatrix(cores, ug)

  wp <- matrix(NA, nrow=length(groups), ncol=length(ug))

  for(i in 1:nrow(wp)){

    use <- unlist(lapply(cores, FUN=function(x){ return( !(i %in% x$core) ) }))
    wp[i,] <- apply(
      matrix(
        unlist(
          lapply(cores, FUN=function(x){
            x$cv.posterior[i,]
          })), ncol=length(ug), byrow=TRUE)[use,] * qm[use,], 2, sum) / apply(qm,2,sum)
  }

  #problem of NaN in class estimation due to NaN in posterior of LDA CV
  ef <- which( apply(wp,1,FUN=function(x){ return(!any(is.nan(x))) }) )
  est <- ug[apply(wp[ef,],1,which.max)]

  return( computeFCR(table(est, groups[ef]), ug)  )

}
