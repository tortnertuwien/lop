extractQualityMatrix <- function(cores, ug){

  qm <- matrix(NA, nrow=length(cores), ncol=length(ug))

  for(i in 1:length(cores)){

    qm[i, ] <- cores[[i]]$quality * cores[[i]]$frequency

  }
  return(qm)

}
