guidedProjections <- function(data, q = 10, osd="OD"){

  getMin <- function(dists, exclude){
    dists[exclude] <- NA
    return(which.min(dists))
  }

  distsLeave1Out <- function(data, fixed, optional, q, which){
    #calculates the distances for each left out observation to
    #the subspace spanned by the remaining observations
    md <- numeric(length(optional))
    for(i in 1:length(optional)){
      md[i] <- getMeasure(ODSD(data, c(fixed, optional[-i]), 1:q), which=which)[optional[i]]
    }
    return(md)
  }

  getStartingObservations <- function(data, number=10){
    dd <- dist(data)
    d11 <- apply(as.matrix(dd), 1, FUN=function(x){ return(sort(x, decreasing=FALSE)[10]) })
    wm <- which.min(d11)[1]
    return(order(as.matrix(dd)[wm,], decreasing=FALSE)[1:number])
  }

  newStartingObsIt <- function(data, start, q, max.iterations=50, which="OD"){

    #similar to newStartingObs; but until new observation is a disavantage

    available <- start
    iterations <- 0

    #first observation:
    ds <- ODSD(data, available, 1:q)
    selection <- getMin(getMeasure(ds, which=which), available)
    d1o <- distsLeave1Out(data, selection, available, q, which=which)


    while( getMeasure(ds, which=which)[selection] < max(d1o) & (iterations<max.iterations)  ){
      #if dist from new obs < than any leave 1 out dist
      # -> remove the largenst leave 1 out obs and replace it by new obs.

      available <- c(available[-which.max(d1o)], selection)

      #calculate new selection and leave 1 out dists
      ds <- ODSD(data, available, 1:q)
      selection <- getMin(getMeasure(ds, which=which), c(available, selection))
      d1o <- distsLeave1Out(data, selection, available, q, which=which)

      iterations <- iterations+1
    }

    #return new starting observations
    return(available)
  }

  ODSD <- function(data, selection, which){

    vf <- apply(data[selection,],2,FUN=function(x){return( (min(x)-max(x))!=0) })
    data <- data[,vf]

    #	for(i in 1:ncol(data)){
    #		crit <- quantile(data[,i],c(0.975,0.025))
    #		data[data[,i]<crit[2],i] <- crit[2]
    #	}

    cen <- apply(data[selection,],2,mean)
    sc <- apply(data[selection,],2,sd)

    data <- scale(data, cen, sc)

    svd <- svd(data[selection,])

    if(length(which)<=length(selection)){
      which.sd <- which[1:(length(selection)-1)]
    }else{
      which.sd <- which
    }

    x.score <- (data%*%svd$v)[,which.sd]

    OD <- sqrt(apply((data-x.score%*%t(svd$v[,which.sd]))^2,1,sum))
    SD <- sqrt(apply( t( t(x.score^2) / ( (svd$d[which.sd])^2 /(length(selection)-1))),1,sum   ))

    return(list(OD=OD,SD=SD))
  }

  getMeasure <- function(odsd, which=c("OD", "SD", "OD*SD")){
    if(which == "SD"){
      comp <- odsd$SD
    }else if(which=="OD"){
      comp <- odsd$OD
    }else if(which =="OD*SD"){
      comp <- odsd$OD * odsd$SD
    }else{
      comp <- odsd$OD * odsd$SD
    }

    return(comp)
    #	return(odsd$OD)
  }







  #phase 0: get first starting observations
  obsInProj <- q
  so <- getStartingObservations(data, q)

  #phase 1: get new - ?better? starting observations; report both lists
  s0 <- so
  selection <- s1 <- newStartingObsIt(data,so,q, which="OD")

  ##remaining/available observations to be selected include the initial starting observations
  available <- (1:nrow(data))[-selection]

  ##initialize the procedure; get first new observation and order the starting observations.
  ds <- ODSD(data,  selection, 1:q)
  m.left <- getMeasure(ds, which=osd)
  no <- getMin(m.left, selection)
  d1o <- distsLeave1Out(data, no, selection, q, which=osd)
  ord <- order(d1o, decreasing=TRUE)
  projections <- rbind(selection[ord], c(selection[ord[-1]],no))
  ds2 <- ODSD(data,  c(selection[ord[-1]],no), 1:q)
  m.right <- getMeasure(ds2, which=osd)
  tsM <- rbind(m.left, m.right)
  tsOD <- rbind(ds$OD, ds2$OD)
  tsSD <- rbind(ds$SD, ds2$SD)

  selection <- c(selection[ord], no)
  available <- (1:nrow(data))[-selection]

  left.right <- numeric()

  #phase 2: add observations to projections - both directions are allowed.
  while(length(available)>1){

    #test whether to add observation at left or right side
    if( min(m.left[available])>min(m.right[available]) ){
      no <- getMin(m.left, selection)
      projObs <- c(no, selection[1:(obsInProj-1)])
      sd <- ODSD(data, projObs, 1:q)
      tsOD <- rbind(sd$OD, tsOD)
      tsSD <- rbind(sd$SD, tsSD)
      m.left <- getMeasure(sd, which=osd)
      tsM <- rbind(m.left, tsM)
      selection <- c(no,selection)
      projections <- rbind(projObs, projections)
      left.right <- c(left.right,1)

    }else{
      no <- getMin(m.right, selection)
      projObs <- c(selection[(length(selection)-obsInProj+2):length(selection)], no)
      sd <- ODSD(data, projObs, 1:q)
      tsOD <- rbind(tsOD,sd$OD)
      tsSD <- rbind(tsSD,sd$SD)
      m.right <- getMeasure(sd, which=osd)
      tsM <- rbind(tsM, m.right)
      selection <- c(selection, no)
      projections <- rbind(projections,projObs)
      left.right <- c(left.right,2)

    }

    available <- (1:nrow(data))[-selection]
    #		print(min(m.left[available])<min(m.right[available]))



  }

  obj <- list(order=selection,  ODs=tsOD, SDs=tsSD, OSD=tsM, projections=projections)
  class(obj) <- "guidedprojections"
  return(obj)


}
