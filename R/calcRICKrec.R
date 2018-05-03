#' Function that calculates recruitment values from Ricker function. For use in \code{\link{calcREF}}
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

calcRICKrec<-function(SRouts,spbioIn)
{
 alpha	<-SRouts[[2]][1]
 beta		<-SRouts[[2]][2]
 Pred		<-spbioIn*exp(alpha-beta*spbioIn)
 return(Pred)
}