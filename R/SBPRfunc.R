#' Calculate spawning-biomass-per-recruit quantities for use in \code{\link{Recruitment}}
#' 
#' @param F Fishing mortality
#' @param Rzero Virgin recruitment
#' @param NatM Natural mortality
#' @param vuln Vector of vulnerability-at-age
#' @param mature Vector of maturity-at-age
#' @param weight Vector of weight-at-age
#' @param MaxAge Maximum age
#'
#' @return List of mature spawning biomass (SBP), spawning biomass per recruit (SBPR), and yield per recruit (YPR).
#'
#' @export

SBPRfunc<-function(F,Rzero,NatM,vuln,mature,weight,MaxAge)
{
  #initial conditions
  N<-matrix(ncol=MaxAge,nrow=50)
  N[1,]<-initialN(Rzero,NatM,MaxAge)
  
  for (j in 2:50)
  {
 for (i in 2:(MaxAge-1))
  N[j,i]<-N[j-1,i-1]*exp(-F*vuln[i-1])*exp(-NatM)

 N[j,MaxAge]<-(N[j-1,MaxAge-1]+N[j-1,MaxAge])*exp(-F*vuln[MaxAge])*exp(-NatM)
 N[j,1]<-Rzero
  }
  SBP<-sum(N[50,]*mature*weight)
  SBPR<-sum(N[50,]*mature*weight)/Rzero
  YPR<-sum(N[50,]*(1-exp(-F*vuln)) *weight)
  list(SBP,SBPR,YPR)
  #return(SBPR)
}
