#' Recruitment function
#' 
#' Applies selected stock-recruitment function using life history parameters. Used in \code{\link{GeMS}}, \code{\link{PlotLifeHistory}}, and \code{\link{ProjPopDym}}.
#' 
#' @param EggsIN Vector of egg mass
#' @param steepnessIN Steepness parameter
#' @param RzeroIN Virgin recruitment
#' @param RecErrIN Recruitment error
#' @param recType "BH" or "Rick"; For Beverton-Holt or Ricker stock-recruitment relationship
#' @param NatMin Natural mortality
#' @param vulnIN Vector of vulnerability-at-age
#' @param matureIN Vector of maturity-at-age
#' @param weightIN Vector of weight-at-age
#' @param LenAtAgeIN Vector of lengths-at-age
#' @param MaxAge Maximum age
#' @param sigmaRin Recruitment variability
#' 
#' @return Vector of recruitment values
#'
#' @export
#'
Recruitment<-function(EggsIN,steepnessIN,RzeroIN,RecErrIN,recType,NatMin,vulnIN,matureIN,weightIN,LenAtAgeIN,MaxAge,sigmaRin)
{
  if(recType=="BH")
  {
 Szero<-unlist(SBPRfunc(0,RzeroIN,NatMin,vulnIN,matureIN,weightIN,MaxAge)[1])	#since there are only a couple possible values, code could be sped up by changing this
 alpha<-(Szero/RzeroIN)*(1-steepnessIN)/(4*steepnessIN)
 beta<-(5*steepnessIN-1)/(4*steepnessIN*RzeroIN)
 Recruit<-(EggsIN/(alpha+beta*EggsIN))*exp(RecErrIN-(sigmaRin^2)/2)
  }
  if(recType=="Rick")
  {
 Szero<-unlist(SBPRfunc(0,RzeroIN,NatMin,vulnIN,matureIN,weightIN,MaxAge)[1])
 beta<-log(5*steepnessIN)/(0.8*Szero)
 alpha<-log(RzeroIN/Szero)+(beta*Szero)
 Recruit<-(EggsIN*exp(alpha-beta*EggsIN))*exp(RecErrIN-(sigmaRin^2)/2)
  }
  return(Recruit)
}