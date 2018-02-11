
bestLocalDiscrimination <- function(data, group, k=NULL){

  if(is.null(k)){
    k.min <- 3/min(table(group))
    k.max <- 0.9
    k <- seq(k.min, k.max, length.out=10)
  }

  dist <- as.matrix(dist(data))
  core <- localDiscrimination(data,group, k=k[1], dist=dist)
  k.opt <- k[1]
  ev <- evaluateCores(core, group)$total
  for(kk in 2:length(k)){
    core.tmp <- localDiscrimination(data,group, k=k[kk], dist=dist)
    ev.tmp <- evaluateCores(core.tmp, group)$total

    if(ev.tmp < ev){
      core <- core.tmp
      ev <- ev.tmp
      k.opt <- k[kk]
    }
  }

  print( paste0("... optimal k set to ", k.opt) )
  return(core)

}
