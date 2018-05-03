#' Function that calculates recruitment values with Beverton-Holt function. For use in \code{\link{calcREF}}
#' Some overlap with \code{\link{BevHolt}}
#'
#' @param SRouts Stock-recruitment parameter outputs
#' @param spbioIn Spawning stock biomass
#' 
#' @return Negative log-likelihood
#'
#' @export
#' @examples
#' \dontrun{
#' 
#' }
#' 
calcBHrec<-function(SRouts,spbioIn)
{
 alpha	<-abs(SRouts[[2]][1])
 beta		<-abs(SRouts[[2]][2])
 Pred		<-alpha*spbioIn/(1+(alpha*spbioIn/beta)) 
 return(Pred)
}