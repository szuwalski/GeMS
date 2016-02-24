#' Beverton-Holt function
#'
#' @param x string of components for BH function
#'
#' @return NLL for fitting BH
#' @export
BevHolt<-function(x)
{
  alpha		<-abs(x[1])
  beta		<-abs(x[2])
  sigma		<-abs(x[3])
  Pred		<-alpha*spbio/(1+(spbio/beta))
  loglike	<-log(1/(sqrt(2*3.141598*sigma^2))) - ((log(Pred)-log(rec))^2)/(2*sigma^2)
  return(sum(-1*loglike))
}