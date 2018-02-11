localDiscrimination <- function(data, group, k, dist=NULL){
  library(MASS)

  if(is.null(dist)){
    dist <- as.matrix(dist(data))
  }

  if(!is.factor(group)){
    stop("group is expected to be a factor variable.")
  }
  if(class(dist)=="dist"){ dist <- as.matrix(dist) }
  if(length(group) != nrow(dist)  ){
    stop("The number of observations differs from the number of labels.")
  }

  cores <- array(list())
  ug <- levels(group)

  if(k<1){
    k0 <- round( k * table(group) ) # changed back to constant k
    k0 <- rep( round(min(k) * table(group) ) , length(ug))
  }else{
    k0 <- rep(k, length(ug))
  }


  for(g in 1:length(ug)){

    gf <- which(group==ug[g])
    d <- dist[gf,gf]
    core.order <- as.factor(apply(d,1, FUN=function(x,lab,k0){
      return( paste(sort(lab[order(x)[1: k0[g] ]]),collapse=",")  )
    }, lab=gf, k0=k0))
    unique.cores <- table(core.order)
    for(i in 1:length(unique.cores)){
      #print(paste0(i, "/", length(unique.cores)))
      core <- as.numeric(unlist(strsplit(names(unique.cores)[i],",")))
      dat <- data[core,]
      loc <- apply(dat,2,mean)
      std <- apply(dat,2,sd)
      vf <- which(std > 0)

      #projection space
      d.s <- scale(data[, vf], center = loc[vf], scale = std[vf])
      svd <- svd(d.s, nu=(length(core)-1), nv = (length(core)-1) )
      d.c <- (d.s %*% svd$v)
      od <- sqrt( apply((d.s - d.c %*% t(svd$v))^2,1, sum) )
      df.core <- scale(as.data.frame(cbind( d.c, od )))

      #quality estimation
      qf  <- (!((1:nrow(df.core)) %in% core)) # ... quality filter, exclude core
      lda <- lda(df.core[qf,], group[qf])

      cv.posterior <- predict(lda, df.core)$posterior
      #cv.posterior <- matrix(NA, nrow=nrow(data), ncol=length(ug))
      #cv.posterior[qf,] <- MASS::lda(df.core[qf,], group[qf], CV=TRUE)$posterior
      q.plus <- q.minus <- numeric( length(ug) )
      for(gg in 1:length(ug)){
        f <- (group == ug[gg])[which(qf)]
        q.plus[ colnames(cv.posterior) == ug[gg] ] <- mean(
          cv.posterior[f, colnames(cv.posterior) == ug[gg] ], na.rm = TRUE )
        q.minus[colnames(cv.posterior) == ug[gg]] <- mean(
          cv.posterior[!f, colnames(cv.posterior) == ug[gg] ], na.rm = TRUE )
      }
      quality <- exp(q.plus - q.minus)
      #fcr <- computeFCR( ug[apply(cv.posterior,1, which.max)], group[qf], levels(group) )

      cores[[length(cores)+1]] <- list(core = core,
                                       group = ug[g],
                                       projection = svd$v,
                                       vf = vf,
                                       location=loc,
                                       scale=std,
                                       core.location=attr(df.core, "scaled:center"),
                                       core.scale=attr(df.core, "scaled:scale"),
                                       frequency = as.numeric(unique.cores[i]),
                                       quality=quality,
                                       cv.posterior=cv.posterior,
                                       lda = lda,
                                       #fcr = fcr
                                       df.core = df.core)
    }
  }

  class(cores) <- "localDiscrimination"
  return(cores)

}
