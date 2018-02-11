weightedLocOut <- function(lo){

  return(
  sapply(1:nrow(lo$ods), FUN=function(x, lo){
    w <- 1/lo$cds[x,]
    weights <- ( w-min(w) ) / sum((w)-min(w) )
    return( sum(weights*lo$ods[x,]) )
  }, lo=lo)
  )

}


weightedLocOut.euc <- function(lo){
  return(
    sapply(1:nrow(lo$ods), FUN=function(x, lo){
      w <- 1/lo$ecds[x,]
      weights <- ( w-min(w) ) / sum((w)-min(w) )
      return( sum(weights*lo$ods[x,]) )
    }, lo=lo)
  )
}


