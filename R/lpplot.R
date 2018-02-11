lpplot <- function(probabilities, class, combination=NA){

  if(max(probabilities) > 1 | min(probabilities) < 0){
    stop("All probabilities must be within 0 and 1.")
  }

  if(length(combination)==1 ){ #show matrix
    par(mar=c(0,0,0,0))
    lpplotmatrix(probabilities, class, labeldist=0.05)

  }else{ #show combination

    lpplotcombination( probabilities, sel=combination,  group = class, angle=90, box=FALSE)
    text(0.05,0.1, combination[1], cex=1)
    text(0.95,0.1, combination[2], cex=1)

  }


}



createRotation <- function(angle){
  return(matrix( c(cos(angle), sin(angle), -sin(angle), cos(angle)),  ncol=2, byrow = TRUE))
}


lpplotcombination <- function(weights, sel, group= NA,  angle=0, box=TRUE, labeldist = 0.05, ...){

  if(length(group)==1){
    group <- as.factor(rep(1,nrow(weights)))
  }

  if(box){
    plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", axes = FALSE, frame.plot=TRUE )
  }else{
    plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", axes = FALSE)
  }

  c <- 0.515
  center <- c(0.5,0.5)
  probs <- c(0.2, 0.4, 0.6, 0.8)

  dir3 <- (c(1,0) %*% createRotation( (angle+0*120)/360*2*pi )) * c
  dir1 <- (c(1,0) %*% createRotation( (angle+1*120)/360*2*pi )) * c
  dir2 <- (c(1,0) %*% createRotation( (angle+2*120)/360*2*pi )) * c
  corner3 <- center + dir3
  corner1 <- center + dir1
  corner2 <- center + dir2
  c12 <- (corner1 + corner2) / 2
  c13 <- (corner1 + corner3) / 2
  c23 <- (corner2 + corner3) / 2

  g13 <- (corner3*(length(levels(group))-2)+corner1)/(length(levels(group)) - 1)
  g23 <- (corner3*(length(levels(group))-2)+corner2)/(length(levels(group)) - 1)
  c3g <- c12 + (corner3-c12)*(length(levels(group))-2)/length(levels(group))

  xc <- c(c13[1], g13[1], c3g[1], g23[1], c23[1], center[1])
  yc <- c(c13[2], g13[2], c3g[2], g23[2], c23[2], center[2])
  pt <- seq(0, 0.5, 0.005)


  polygon(xc, yc, col="#e0e0e0", border=TRUE, lty=5)

  points( rbind(corner1, corner2), type='l')
  points( rbind(c13, corner3), type='l', lty=1)
  points( rbind(corner1, c13), type='l')
  points( rbind(c23, corner3), type='l', lty=1)
  points( rbind(corner2, c23), type='l')

  u <- levels(group)
  w <- weights/apply(weights,1,sum)
  if(ncol(w)==3){
    df <- data.frame( w[,sel[1]], w[,sel[2]], w[,-sel])
  }else{
    df <- data.frame( w[,sel[1]], w[,sel[2]], apply(w[,-sel],1,sum))
  }
  names(df) <- c(u[sel], "r")

  #grid
  colp05 <- gray(0.5)


  p12c <- (corner1+corner2)/2
  center <- (corner1+corner2+corner3)/3
  pcp <- (p12c - (corner3-p12c))
  for(pp in probs){
    p13p <- (pp*corner1+(1-pp)*corner3)
    p23p <- (pp*corner2+(1-pp)*corner3)
    #if(pp>0.5){
    if(TRUE){
      pc1 <- (pp*corner1+(1-pp)*corner2)
      pc2 <- (pp*corner2+(1-pp)*corner1)
    }else{
      pc1 <- pc2 <- (pp*pcp+(1-pp)*corner3)
    }
    points( rbind(p23p, pc2), type='l', lty=2, col=colp05)
    points( rbind(p13p, pc1), type='l', lty=2, col=colp05)
    dir <- (corner2-corner1); dir <- dir/sum(abs(dir))
    text( p13p - dir * labeldist, as.character(pp), col=colp05, cex=0.75, srt=angle-90)
    text( p23p + dir * labeldist, as.character(pp), col=colp05, cex=0.75, srt=angle-90)

  }
  #p12c
  points( rbind(c3g, c12), type='l', lty=5)
  #points( rbind(c12, center), type='l', lty=5)
  #points( rbind(c13, center), type='l', lty=5)
  #points( rbind(c23, center), type='l', lty=5)


  transf <- function(wx,dir1,dir2,dir3, center){
    return(center + wx[1]*dir1 + wx[2]*dir2 + wx[3]*dir3)
  }
  td <- t(apply(df,1,FUN=transf, dir1=dir1, dir2=dir2, dir3=dir3, center=center))

  pch <- rep(1,length(group))
  pch[group == u[sel[1]]] <- 3
  pch[group == u[sel[2]]] <- 4
  points(td, pch=pch, col=as.numeric(group))

}


lpplotmatrix <- function(weights, group, labeldist=0.05, angle=45, ...){

  u <- levels(group)
  d <- length(u)
  nr <- (d-1)
  ratio <- 5
  c <- 0.5

  layout <- matrix(NA, nrow = nr*ratio+1, ncol=nr*ratio+1)
  cnt <- 1
  for(i in 1:nr){
    for(j in 1:(nr+1-i)){
      cstart <- (i-1)*ratio+2 + (j-1)*ratio
      cend <- (i)*ratio+1 + (j-1)*ratio
      rstart <- (i-1)*ratio+1
      rend <- (i)*ratio
      layout[rstart:rend, cstart:cend] <- cnt
      cnt <- cnt+1
    }
  }

  layout[1:ratio, 1] <- cnt; cnt <- cnt+1
  for(i in 2:nr){
    rstart <- (i-1)*ratio+1
    rend <- (i)*ratio
    cstart <- (i-2)*ratio+2
    cend <- (i-1)*ratio+1
    layout[rstart:rend, cstart:cend] <- cnt
    cnt <- cnt+1
  }
  layout[nrow(layout), (ncol(layout)-ratio+1):ncol(layout)] <- cnt; cnt <- cnt+1


  # for(i in 1:nr){
  #   rstart <- (i-1)*ratio+2
  #   rend <- (i)*ratio+1
  #   layout[ncol(layout), rstart:rend] <- cnt
  #   cnt <- cnt+1
  # }
  # layout[nrow(layout), 1] <- 0
  cstart <- 2
  cend <- (nr-2)*ratio+1
  rstart <- (nr-1)*ratio+1
  rend <- (nr)*ratio
  layout[rstart:rend, cstart:cend] <- cnt

  layout[is.na(layout)] <- 0
  layout(layout)

  for(i in 1:d){ for(j in 1:d){
    if(j > i){
      lpplotcombination(weights, group, sel=c(i, j), angle=45, labeldist=labeldist)
    }
  }}
  # createRotation <- function(angle){
  #   return(matrix( c(cos(angle), -sin(angle), sin(angle), cos(angle)),  ncol=2, byrow = TRUE))
  # }

  dir1 <- c(0.5, 0.5) + (c(1,0) %*% createRotation( (angle+1*120)/360*2*pi )) * c
  dir2 <- c(0.5, 0.5) + (c(1,0) %*% createRotation( (angle+2*120)/360*2*pi )) * c


  plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", axes = FALSE)
  text(0.5, 2/3, u[1], cex=2)

  for(i in 2:(d-1)){
    plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", axes = FALSE)
    text(2/3, 2/3, u[i], cex=2)
  }
  plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", axes = FALSE)
  text(2/3, 0.5, u[length(u)], cex=2)


  plot(NA, xlim=c(0,1), ylim=c(0,1), type="n", xlab="", ylab="", axes = FALSE)
  legend("bottomleft", pch=c(rep(1, d), 3,4),
         col=c(u,grey(0.25), grey(0.25)),
         legend=c(paste0("Class ", u), "1st selected class", "2nd selected class"),
         cex=1.5, bty = "n")

}



