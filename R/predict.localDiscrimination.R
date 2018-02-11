predict.localDiscrimination <- function(cores, data, ...){

  ug <- levels(as.factor(unlist(lapply(cores, FUN=function(x){return(x$group)}))) )


  posteriors <- array(NA, c(length(cores), nrow(data),length(ug)))
  qm <- extractQualityMatrix(cores, ug)
  freq <- unlist(lapply(cores, FUN=function(x){ return(x$frequency) } ))

  for(i in 1:length(cores)){

    #projection space
    d.s <- scale(data[, cores[[i]]$vf],
                 center = cores[[i]]$location[cores[[i]]$vf],
                 scale = cores[[i]]$scale[cores[[i]]$vf])
    d.c <- (d.s %*% cores[[i]]$projection)
    od <- sqrt( apply((d.s - d.c %*% t(cores[[i]]$projection))^2,1, sum) )
    df.core <- scale(as.data.frame(cbind( d.c, od )),
                     center = cores[[i]]$core.location,
                     scale= cores[[i]]$core.scale)
    posteriors[i,,] <- predict(cores[[i]]$lda, df.core)$posterior
  }

  wp <- wp.raw <- matrix(NA, nrow=nrow(data), ncol=length(ug))
  for(i in 1:nrow(wp)){
    wp[i,] <- apply(posteriors[,i,] * qm,2,sum) / apply(qm,2,sum)
    wp.raw[i,] <- apply(posteriors[,i,],2,mean)
  }

  return(
    list(weights=wp,
         class=ug[apply(wp,1,which.max)],
         weights.raw = wp.raw,
         class.raw = ug[apply(wp.raw,1,which.max)])
  )

}
