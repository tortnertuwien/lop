
 getLocalDescription <- function(data, center, dist=NULL, k=20, alpha=0.5){

  if(is.null(dist)){
    dist <- as.matrix(dist(data))
  }else{
    dist <- as.matrix(dist)
  }

  knn <- order( dist[center,], decreasing = FALSE )[2:(k+1)]

  sortPossibeObs <- apply(dist[knn,knn],1,
                          FUN=function(x, alpha){
                            sort(x, decreasing=FALSE)[ceiling(k*alpha)]
                          }, alpha=alpha)
  selcore <- knn[which.min(sortPossibeObs)]
  selection <- order( dist[selcore,], decreasing = FALSE )[1:ceiling(k*alpha)]

  loc <- apply(data[selection,],2,mean)
  std <- apply(data[selection,],2,sd)
  vf <- which(std>0)
  data <- scale(data, loc, std)[,vf]

  svd <- svd(data[selection,])

  dim.proj <- min((ceiling(k*alpha)-1), ncol(data))
  which.sd <- 1:dim.proj

  x.score <- (data%*%svd$v)[,which.sd]

  if(dim.proj < ncol(data)){
    OD <- sqrt(apply((data-x.score%*%t(svd$v[,which.sd]))^2,1,sum))
  }else{
    OD <- OD.k <- rep(0,nrow(data))
  }
  CD <- sqrt(apply( t( t(x.score^2) / ( (svd$d[which.sd])^2 /(length(selection)-1))),1,sum   ))
  ECD <- sqrt( apply(x.score^2,1,sum) )

  return(list(od=OD, cd=CD,
              ecd=ECD,
              core=selection, knn=knn, center=center))
 }



