#' StrLine
#'
#' @param x
#'
#' @return no idea
#' @export
StrLine<-function(x)
{
  Pred		<-x[1]
  sigma	<-x[2]
  loglike	<- log(1/(sqrt(2*3.141598*sigma^2))) - ((log(Pred)-log(rec))^2)/(2*sigma^2)
  return(sum(-1*loglike))
}