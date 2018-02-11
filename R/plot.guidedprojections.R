plot.guidedprojections <- function(gp, col=NA, lty=NA, ...){

  matrix <- t(gp$OSD)
  xl <- c(0,ncol(matrix))
  yl <- c(0,max(matrix))

  xlab="Projection Sequence"
  ylab="OSD"

  #par(mar=c(4.5,4.5,2,1), mgp = c(3, 1, 0))
  plot(NA, xlab=xlab, ylab=ylab, xlim = xl, ylim=yl, bty="n",axes=FALSE)
  axis(1, col = "white", tcl = 0)
  axis(2, col = "white", tcl = 0)
  grid(lty="solid")

  for(i in 1:nrow(matrix)){
    points(1:ncol(matrix), matrix[i,], type='l', col=col[i], lty=lty[i])
  }

}




plot.guidedprojections2 <- function(gp, col=NA, coll=NA){
  library(ggplot2)
  df <- as.data.frame(gp$OSD)
  colnames(df) <- 1:ncol(df)

  #coll <- c("#00ff00","#ff0000")
  if(is.na(col[1])){
    col <- rep(1,ncol(gp$OSD))
  }else if(length(col) != ncol(gp$OSD)){
    col <- rep(1,ncol(gp$OSD))
    warning("The number of observations does not match the number of provided colors.")
  }
  dat <- data.frame(x=NA, y=NA, grp=NA, col=NA)
  for(i in 1:ncol(gp$OSD)){
    dat <- rbind(dat, data.frame(x=1:nrow(gp$OSD), y=df[,i], grp=rep(i,nrow(gp$OSD)), col=rep(col[i], nrow(gp$OSD))))
  }
  dat <- dat[-1,]

  gg<-ggplot(dat, aes(x=x, y=y, group=grp))+
    xlab("Projection Sequence") +
    ylab("OSD") +
    theme(axis.line=element_blank(),
          axis.ticks=element_blank(),
          legend.position="none",
          panel.grid.major = element_line(colour = "#cccccc"),
          panel.grid.minor = element_line(colour = "#eeeeee"),
          panel.background = element_rect(fill = "white")
    )+geom_line(colour=dat$col)
  plot(gg)

}
