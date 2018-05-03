#' Negative log-likelihood
#' 
#' Used in \code{\link{StockRecCheck}}.
#' 
#' @param x Vector of length two, c(predicted value and sigma)
#' @param rec Vector of recruitment values
#' 
#' @return Negative log-likelihood value
#'
StrLine<-function(x,rec)
{
 Pred		<-x[1]
 sigma	<-x[2]
 loglike	<- log(1/(sqrt(2*pi*sigma^2))) - ((log(Pred)-log(rec))^2)/(2*sigma^2) 
 return(sum(-1*loglike))
}