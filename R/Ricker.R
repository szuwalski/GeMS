#' Ricker function
#' 
#' @return Negative log-liklihood from fitting to Ricker stock-recruitment model
#'
#' @export
#' 
Ricker<-function(x)
{
alpha		<-x[1]
beta		<-x[2]
sigma		<-abs(x[3])
Pred		<-spbio*exp(alpha-beta*spbio)
loglike	<-log(1/(sqrt(2*3.141598*sigma^2))) - ((log(Pred)-log(rec))^2)/(2*sigma^2) 
return(sum(-1*loglike))
}