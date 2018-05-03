#' Function that calculates negative log-likelihood for Beverton-Holt recruitment function. For use in \code{\link{StockRecCheck}}
#' Some overlap with \code{\link{calcBHrec}}
#'
#' @param x Vector of Beverton-Holt parameters c(alpha, beta, sigma)
#' @param spbio Spawning Stock Biomass
#' 
#' @return Negative log-likelihood
#'
#' @export
#' @examples
#' \dontrun{
#' x <- c(alpha, beta, 0.001)
#' BevHolt(x, 100000)
#' }

BevHolt<-function(x,spbio)
{
alpha		<-abs(x[1])
beta		<-abs(x[2])
sigma		<-abs(x[3])
Pred		<-alpha*spbio/(1+(spbio/beta))
loglike	<-log(1/(sqrt(2*3.141598*sigma^2))) - ((log(Pred)-log(rec))^2)/(2*sigma^2) 
return(sum(-1*loglike))
}