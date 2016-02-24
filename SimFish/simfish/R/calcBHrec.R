#' calculate recruits per beverton-holt
#'
#' @param SRouts spawner recruit parameters
#' @param spbioIn spawning biomass
#'
#' @return predicted recruits
#' @export
calcBHrec<-function(SRouts,spbioIn)
{
  alpha	<-abs(SRouts[[2]][1])
  beta		<-abs(SRouts[[2]][2])
  Pred		<-alpha*spbioIn/(1+(alpha*spbioIn/beta))
  return(Pred)
}
