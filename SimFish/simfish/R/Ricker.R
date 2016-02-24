#' Fit ricker model
#'
#' @param x string of parameters for ricker model
#'
#' @return nll for ricker
#' @export
Ricker<-function(x)
{
  alpha		<-x[1]
  beta		<-x[2]
  sigma		<-abs(x[3])
  Pred		<-spbio*exp(alpha-beta*spbio)
  loglike	<-log(1/(sqrt(2*3.141598*sigma^2))) - ((log(Pred)-log(rec))^2)/(2*sigma^2)
  return(sum(-1*loglike))
}