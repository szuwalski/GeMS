#' Calculate recruitment
#'
#' \code{Recruitment} returns the recruits for a given functional form and SSB
#' @param EggsIN eggs
#' @param steepnessIN steepness
#' @param RzeroIN eq recruitment
#' @param RecErrIN recruitment error
#' @param recType form of recruitment
#' @param NatMin natural mortality
#' @param vulnIN fishing vulnerability
#' @param matureIN maturity at age
#' @param weightIN weight at age
#' @param LenAtAgeIN length at age
#' @param MaxAge maximim age
#'
#' @return returns the recruits
#' @export
Recruitment<-function(EggsIN,steepnessIN,RzeroIN,RecErrIN,recType,NatMin,vulnIN,matureIN,weightIN,LenAtAgeIN, MaxAge)
{
  if(recType=="BH")
  {
    Szero<-unlist(SBPRfunc(F = 0, Rzero = RzeroIN,NatM = NatMin,vuln = vulnIN,
                           mature = matureIN, weight= weightIN, LenAtAge = LenAtAgeIN,
                           MaxAge = MaxAge)[1])	#since there are only a couple possible values, code could be sped up by changing this
    alpha<-(Szero/RzeroIN)*(1-steepnessIN)/(4*steepnessIN)
    beta<-(5*steepnessIN-1)/(4*steepnessIN*RzeroIN)
    Recruit<-(EggsIN/(alpha+beta*EggsIN))*exp(RecErrIN)
  }
  if(recType=="Rick")
  {
    Szero<-unlist(SBPRfunc(0,RzeroIN,NatMin,vulnIN,matureIN,weightIN, LenAtAgeIN,MaxAge)[1])
    beta<-log(5*steepnessIN)/(0.8*Szero)
    alpha<-log(RzeroIN/Szero)+(beta*Szero)
    Recruit<-(EggsIN*exp(alpha-beta*EggsIN))*exp(RecErrIN)
  }
  return(Recruit)
}