#' Calculate initial biomass using \eqn{Recruitment * exp(-NaturalMortality)}
#' 
#' Used in \code{\link{ProjPopDym}}, \code{\link{FindF35}}, \code{\link{SBPRfunc}}
#' 
#' @param Rzero Virgin recruitment
#' @param NatM Natural mortality
#' @param inAge Maximum age
#'
#' @return Vector of initial biomass
#'
#' @export
#'
initialN<-function(Rzero,NatM,inAge)
{
  test		<-rep(0,100)
  test[1]	<-Rzero[1]
  for (i in 2:100)
   {test[i]	<-test[i-1]*exp(-NatM)}
  virInit2	<-c(test[1:(inAge-1)],sum(test[inAge:100]))
  return(virInit2)
}