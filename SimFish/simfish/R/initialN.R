#' Calculate initial numbers
#'
#' @param Rzero recruits at EQ
#' @param NatM natural mortality
#' @param inAge no idea
#'
#' @return initial numbers
#' @export
initialN<-function(Rzero,NatM,inAge)
{
  test		<-rep(0,100)
  test[1]	<-Rzero[1]
  for (i in 2:100)
  {test[i]	<-test[i-1]*exp(-NatM)}
  virInit2	<-c(test[1:(inAge-1)],sum(test[inAge:100]))
  return(virInit2)
}