#' Calculate ricker recruits
#'
#' @param SRouts parameters of ricker recruitment
#' @param spbioIn spawning biomass
#'
#' @return predicted ricker recruits
#' @export
calcRICKrec<-function(SRouts,spbioIn)
{
  alpha	<-SRouts[[2]][1]
  beta		<-SRouts[[2]][2]
  Pred		<-spbioIn*exp(alpha-beta*spbioIn)
  return(Pred)
}
