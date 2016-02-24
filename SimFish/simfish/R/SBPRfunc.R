#' Calc Spawning biomass per recruit
#'
#' \code{SBPRfunc} calculate spaning biomass per recruit based on
#' parameters
#' @param Rzero recruitment at eq
#' @param NatM natural mortality
#' @param vuln vulnerability to fishing
#' @param mature maturity at age
#' @param weight weight at age
#' @param LenAtAge length at age

SBPRfunc<-function(F,Rzero,NatM,vuln,mature,weight,LenAtAge, MaxAge)
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
